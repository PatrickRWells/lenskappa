from lenskappa.analysis.transformation import Transformation
from lenskappa.utils.attach_ms_wlm import attach_wlm
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


def delegate_weights(weights_list, nweights):
    combs = list(combinations(weights_list, nweights))
    return combs

class load_wnc(Transformation):
    def __call__(self, wnc_path: str):
        path = Path(wnc_path)
        if not path.exists():
            raise FileNotFoundError(f"No file found at {str(path)}")
        print("reading file")
        return pd.read_csv(path)
    
class build_wnc_distribution(Transformation):

    def __call__(self, wnc: pd.DataFrame, weights = ["gal", "massoverr", "zoverr"], bins_per_dim = 100, weights_min = 0, weights_max = 3, *args, **kwargs):
        selected_weights = wnc[weights].to_numpy()
        if type(weights_min) == dict:
            bounds = [(weights_min[w], weights_max[w]) for w in weights]
        else:
            bounds = [(weights_min, weights_max) for _ in weights]
        hist, edges = np.histogramdd(selected_weights, bins_per_dim, bounds, density=True)
        return hist, edges
    
class load_ms_wnc(Transformation):
    def __call__(self, ms_wnc_path):
        print(f"Loading and normalize MS weights")
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
    def __call__(self, ms_wnc, z_s, wlm_path, ms_wnc_path, threads=1):
        all = []
        missing = {}
        for field, weights in ms_wnc.items():
            if 'gamma' in weights.columns and 'kappa' in weights.columns:
                all.append(weights)
            else:
                missing.update({field: weights})
        if missing:
            remaining = attach_wlm(missing, z_s, Path(wlm_path), Path(ms_wnc_path), threads)
            all = all + [remaining]
        all_weights = pd.concat(all)
        return all_weights
    
class partition_ms_weights(Transformation):
    def __call__(self, ms_weights_wwlm, weights, wnc_distribution):
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
    def __call__(self, ms_weights_wwlm, ms_weight_partitions, wnc_distribution, kappa_bins = None, output_path = None, threads=1):
        if kappa_bins is None:
            kappa_bins = np.linspace(-0.2, 0.4, 1000)
        pdf = np.zeros_like(kappa_bins[:-1])
        weight_pdf = wnc_distribution[0]
        kappas = ms_weights_wwlm["kappa"].to_numpy()
        idxs = [np.array(idx) for idx in np.ndindex(weight_pdf.shape)]
        f_ = partial(compute_single_pdf, ms_partitions = ms_weight_partitions, kappas = kappas, kappa_bins = kappa_bins)
        with mp.Pool(threads) as p:
            chunksize = len(idxs) // 100
            results = list(tqdm.tqdm(p.imap(f_, idxs, chunksize=chunksize), total=len(idxs)))
            
        results = np.array([r*weight_pdf[tuple(idxs[i])] for i, r in enumerate(results)])
        pdf = np.sum(results, axis=0)
        if output_path is not None:
            output = {"bins": kappa_bins, "pdf": pdf}
            with open(output_path, "wb") as f:
                pickle.dump(output, f)

        return kappa_bins, pdf
    
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
