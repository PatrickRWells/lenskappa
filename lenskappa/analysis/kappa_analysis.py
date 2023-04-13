from lenskappa.analysis.transformation import Transformation
from lenskappa.utils.attach_ms_wlm import attach_wlm, get_redshift_plane
import pandas as pd
from pathlib import Path
import re
import numpy as np
import multiprocessing as mp
from functools import partial
import numba
import tqdm
from itertools import combinations
import pickle
import logging
from dask.distributed import get_client, secede, rejoin
from dask import array
import time
import psutil
import os
import random
"""
This analysis represents a single kappa inference,
using a particular set of weights.

Parameters
----------

wnc_path: Path
    The path to the weighted number counts for the lens system. Should
    be stored as a CSV
ms_wnc_path: Path
    The path to the weighted number counts from the simulation. Should
    be a folder with (potentially) several files
weights: list
    The list of weights to use in the inference
z_s: float
    The source redshift. 
weights_min: float or dict, default = 0
    The lowest value of weights to include in the inference. If a float,
    this value will be used for all weights. If dictionary, should
    have an entry for each weight and look like {"weight-name": min-value}
weights_max: float or dict, default = 3.0
    The highest value of the weights to include in the inference, using the same
    format as weights_min
bins_per_dim: int, default = 100
    The number of bins to partition each weight into. The actual number
    of bins will be bins_per_dim^n_weights
kappa_bins: list[float], defualt = np.linspace(-0.2, 0.4, 1000)
    The bin edges for the final kappa histogram. 



"""

def delegate_weights(weights_list, nweights):
    combs = list(combinations(weights_list, nweights))
    return combs

def setup(config):
    ret = {}
    weight_path = config["parameters"]["wnc_path"]
    spec = weight_path.stem
    spec = spec.split("_")
    spec = {"limmag": spec[0], "aperture": spec[1]}
    if "redshift_plane" not in config["parameters"]:
        redshift_plane = get_redshift_plane(config["parameters"]["z_s"])
        ret.update({"redshift_plane": redshift_plane})
    ret.update({"spec": spec})
    return ret



class load_wnc(Transformation):
    def __call__(self, wnc_path: str, logger_: logging.Logger, **kwargs):
        path = Path(wnc_path)
        if not path.exists():
            raise FileNotFoundError(f"No file found at {str(path)}")
        return pd.read_csv(path)
    
class build_wnc_distribution(Transformation):

    def __call__(self, wnc: pd.DataFrame, spec: dict, logger_: logging.Logger, name, weights = ["gal", "massoverr", "zoverr"], bins_per_dim = 100, weights_min = 0, weights_max = 3, *args, **kwargs):
        weight_str = ", ".join(weights)
        logger_.info(f"SYSTEM {name}: LIMMAG {spec['limmag']}: APERTURE {spec['aperture']}: Building distribution for weighted number counts for weight combination {weight_str}" )
        selected_weights = wnc[weights].to_numpy()
        if type(weights_min) == dict:
            bounds = [(weights_min[w], weights_max[w]) for w in weights]
        else:
            bounds = [(weights_min, weights_max) for _ in weights]
        hist, edges = np.histogramdd(selected_weights, bins_per_dim, bounds, density=True)
        return hist, edges
    
class load_ms_wnc(Transformation):
    def __call__(self, ms_wnc_path, logger_: logging.Logger, **kwargs):
        #logger_.info(f"Loading and normalize MS weights")
        path = Path(ms_wnc_path)
        file_paths = [f for f in path.glob("*.csv")]
        data = {}
        for path in file_paths:
            name = path.name
            field = re.findall(r"\d", name)
            key = tuple(int(f) for f in field)
            field_data = pd.read_csv(path)
            data.update({key: self.normalize_ms_weights(field_data)})
        return data


    def normalize_ms_weights(self, weights):
        columns_to_skip = ["ra", "dec", "kappa", "gamma"]
        output_weights = weights.replace([-np.inf, np.inf], np.nan)
        output_weights = output_weights.dropna()
        for col in output_weights.columns:
            if col in columns_to_skip:
                continue
            output_weights[col] = weights[col] / weights[col].median()
        return output_weights

class attach_ms_wlm(Transformation):
    def __call__(self, ms_wnc, z_s, wlm_path, ms_wnc_path, logger_: logging.Logger, redshift_plane = None, **kwargs):
        all = []
        missing = {}
        for field, weights in ms_wnc.items():
            if 'gamma' in weights.columns and 'kappa' in weights.columns:
                all.append(weights)
            else:
                missing.update({field: weights})
        if missing:
            remaining = attach_wlm(missing, z_s, redshift_plane, Path(wlm_path), Path(ms_wnc_path))
            all = all + [remaining]
        all_weights = pd.concat(all)
        return all_weights
    
class partition_ms_weights(Transformation):
    def __call__(self, ms_weights_wwlm, weights, wnc_distribution, logger_: logging.Logger, **kwargs):
        edges = wnc_distribution[1]
        ms_weight_values = ms_weights_wwlm[weights].to_numpy()
        indices = np.empty(ms_weight_values.shape[1], dtype=object)
        for i, column in enumerate(ms_weight_values.T):
            #finds which box each los belongs in
            ixs = np.digitize(column, edges[i]) - 1
            #digitize gives you the index such that edges[i-1] <= val <= edges[i]
            #so we substract one so the index matches up with the bin centers
            #defined in kappa_pdf
            indices[i] = ixs

    #this gives us a 2d array, where each row is the set of indices
    #that correctly bins the given sample in the millennium simulation
        return np.array(list(zip(*indices)), dtype=int)


class compute_pdfs(Transformation):
    def __call__(self, ms_weights_wwlm, spec, ms_weight_partitions, wnc_distribution, weights, logger_: logging.Logger, kappa_bins = None, output_path = None, name=None):
        if kappa_bins is None:
            kappa_bins = np.linspace(-0.2, 0.4, 1000)
        pdf = np.zeros_like(kappa_bins[:-1])
        weight_pdf = wnc_distribution[0]
        kappas = ms_weights_wwlm["kappa"].to_numpy()
        weights_str = ", ".join(weights)
        idxs = [np.array(idx) for idx in np.ndindex(weight_pdf.shape)]
        #f_ = partial(compute_single_pdf, ms_partitions = ms_weight_partitions, kappas = kappas, kappa_bins = kappa_bins)
        client = get_client()
        nworkers = len(client.scheduler_info()["workers"])
        nper = len(idxs) // nworkers
        random.shuffle(idxs)
        results = []
        for i in range(nworkers):
            if i == (nworkers - 1):
                idx_range = idxs[nper*i: -1]
            else:
                idx_range = idxs[nper*i: nper*(i+1)]
            f_ = partial(compute_pdf_range, idx_list = idx_range, weight_pdf = weight_pdf, ms_partitions = ms_weight_partitions, kappas = kappas, kappa_bins = kappa_bins)
            results.append(client.submit(f_))
        secede()
        logger_.info(f"SYSTEM {name}: LIMMAG {spec['limmag']}: APERTURE {spec['aperture']}: computing kappa histogram for weight combination {weights_str}")
        results = np.array(client.gather(results))
        rejoin()
        pdf = np.sum(results, axis=0)
        if output_path is not None:
            output = {"bins": kappa_bins, "pdf": pdf}
            with open(output_path, "wb") as f:
                pickle.dump(output, f)

        return kappa_bins, pdf


def compute_pdf_range(idx_list, weight_pdf, ms_partitions, kappas, kappa_bins):
    results = np.array([compute_single_pdf(idx, ms_partitions, kappas, kappa_bins) for idx in idx_list])
    results = np.array([result*weight_pdf[tuple(idx_list[i])] for i, result in enumerate(results)])
    results = np.sum(results, axis=0)
    return results

@numba.njit
def compute_single_pdf(idx, ms_partitions, kappas, kappa_bins):
        mask = np.ones(len(ms_partitions), dtype=np.bool8)
        bi_mask = (ms_partitions == idx)
        for i in range(bi_mask.shape[1]):
            mask = np.logical_and(mask, bi_mask[:,i])
        if np.any(mask):
            histogram, _ = np.histogram(kappas[mask], bins = kappa_bins)
            histogram = histogram/np.count_nonzero(mask)
        else:
            histogram = np.zeros(len(kappa_bins) - 1, dtype=np.float64)
        return histogram
