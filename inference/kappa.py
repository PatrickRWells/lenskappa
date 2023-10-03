from ast import arg
from concurrent.futures import process
import itertools
import logging
import pathlib
from statistics import median
from typing import Tuple, Type
import numpy as np
import pandas as pd
import re
import multiprocessing as mp
from lenskappa.utils import attach_wlm
import astropy.units as u
import random
from scipy import stats
from tqdm import tqdm
from functools import partial
import numba
from copy import copy
import tqdm

dist = [float, float, float]


def kappa_pdf(
        lens_wnc: pathlib.Path, ms_wnc: pathlib.Path, z_s: float, 
        weights = ["gal", "massoverr", "zoverr"], min_kappa = -0.2,
        max_kappa = 0.4, kappa_bins = 1000, *args, **kwargs
        ):
    ms_weights = load_ms_weights(ms_wnc)
    if "kappa" not in ms_weights.columns:
        ms_weights = attach_wlm(ms_wnc, z_s)
    
    field_weights = pd.read_csv(lens_wnc)
    field_weights = field_weights.replace([-np.inf, np.inf], np.nan)
    field_weights = field_weights.dropna()
    cols_to_skip = ["Unnamed: 0", "ra", "dec", "kappa", "gamma"]
    normalized_ms_weights = ms_weights.copy()
    
    for col in ms_weights.columns:
        if col in cols_to_skip:
            continue
        normalized_ms_weights[col] = normalized_ms_weights[col]/ms_weights[col].median()
        
    print("Building likelihood distribution from weighted number counts...")
    weight_pdf, edges = build_distribution(field_weights, weights=weights, *args, **kwargs)
    bin_centers = np.empty(len(edges), dtype = object)
    for i, ed in enumerate(edges):
        bin_centers[i] = (ed[:-1] + ed[1:])/2
    ms_weight_values = normalized_ms_weights[weights].to_numpy(dtype=np.float32)
    kappa = normalized_ms_weights["kappa"].to_numpy(dtype=np.float32)
    print("Partitioning samples from the Millennium Simulation")
    ms_partitions = partition_ms(ms_weight_values, edges)


    kappa_bins = np.linspace(min_kappa, max_kappa, kappa_bins, dtype=np.float32)
    print("Computing histograms.... This may take some time")
    print("If the number of weights considered >= 3, this bar may only update every minute or two.")
    pdf = compute_pdfs(weight_pdf, ms_partitions, kappa, kappa_bins, *args, **kwargs)

    mean_bins = (kappa_bins[:-1] + kappa_bins[1:]) / 2
    return mean_bins, pdf

def compute_pdfs(weight_pdf, ms_partitions, kappas, kappa_bins, threads = 1):
    pdf = np.zeros_like(kappa_bins[:-1])
    idxs = [np.array(idx) for idx in np.ndindex(weight_pdf.shape)]
    f_ = partial(compute_single_pdf, ms_partitions = ms_partitions, kappas = kappas, kappa_bins = kappa_bins)
    with mp.Pool(threads) as p:
        chunksize = len(idxs) // 100
        results = list(tqdm.tqdm(p.imap(f_, idxs, chunksize=chunksize), total=len(idxs)))
        
    results = np.array([r*weight_pdf[tuple(idxs[i])] for i, r in enumerate(results)])
    pdf = np.sum(results, axis=0)

    return pdf

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


def load_ms_weights(path: pathlib.Path):
    files = [f for f in path.glob('*.csv') if not f.name.startswith('.')]
    if len(files) != 64:
        print(f"Warning: Expected 64 MS weight files but only found {len(files)}")
    dfs = [pd.read_csv(f) for f in files]
    return pd.concat(dfs)

def build_distribution(field_weights: pd.DataFrame, weights = ["gal", "massoverr", "zoverr"], bins_per_dim = 100, lower_bound = 0, upper_bound = 3, *args, **kwargs):
    selected_weights = field_weights[weights].to_numpy()
    hist, edges = np.histogramdd(selected_weights, bins_per_dim, [(lower_bound, upper_bound) for _ in weights], density=True)
    return hist, edges

def partition_ms(ms_weight_values: np.array, edges):
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

