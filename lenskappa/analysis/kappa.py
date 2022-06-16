from ast import arg
from concurrent.futures import process
import itertools
import logging
import pathlib
from statistics import median
from typing import Tuple, Type
from lenskappa.analysis.kappa_old import compute_single_histogram
from lenskappa.utils import SurveyDataManager
import numpy as np
import pandas as pd
import re
import multiprocessing as mp
from lenskappa.datasets.surveys.ms.ms import millenium_simulation
import astropy.units as u
import random
from scipy import stats
from tqdm import tqdm
from functools import partial

DATA_MANAGER = SurveyDataManager("ms")
dist = [float, float, float]

def kappa_pdf(dists: dist, simulation: str = "ms", *args, **kwargs):
    logger = logging.getLogger("kappa_inference")
    if simulation == "ms":
        config, data = setup_ms(logger, *args, **kwargs)
    
    config['threads'] = kwargs.get("threads", 1)
    dists = normalize_all(dists, data)
    bins = get_bins(dists, data, *args, **kwargs)
    combos = get_bin_combos(bins, *args, **kwargs)
    random.shuffle(combos) #To evenly distribute work
    l = len(combos)
    workers = config['threads'] - 1
    nperthread = round(len(combos)/workers)
    thread_bins = [i*nperthread for i in range(workers)]
    thread_bins.append(l)
    binned_combos = [combos[b: thread_bins[i+1]] for i, b in enumerate(thread_bins[:-1])]
    processes = []
    #mgr = mp.Manager()
    #ns = mgr.Namespace()
    queue = mp.Queue()
    #ns.bins = bins
    #ns.data = data
    print("Delegating...")
    #compute_histogram_range(combos, bins, data, None)
    for index, b in enumerate(binned_combos):
        print(f"Thread {index} gets {len(b)} histograms")
        p = mp.Process(target=compute_histogram_range, args=(b, bins, data, queue))
        p.start()
        processes.append(p)
    hists = []
    for p in processes:
        hists.append(queue.get())

    for p in processes:
        p.join()
    hists = list(itertools.chain(*hists))
    hists_values = [h[0] for h in hists]
    bins = hists[0][1]
    return bins, np.sum(hists_values, axis = 0)
def compute_histogram_range(combos, bins, data, queue):
    vals = []
    for i,c in enumerate(combos):
        if i%100 == 0:
            print(f"Completed {i}")
            print(f"Found {len(vals)} bins with LOS in them")
        values = compute_single_histogram(c, bins, data)
        if values:
            vals.append(values)
    print("Putting in queue!")
    queue.put(vals)

def compute_single_histogram(combo, bins, data):
    gauss = stats.norm(0, 1)
    gauss_factor = 1.0
    weights = data['ms_weights']
    names = list(bins.keys())
    mask = np.ones(len(weights), dtype=bool)
    for index, name in enumerate(names):
        d = bins[name][combo[index]]
        mask = mask&d['mask']
        gauss_factor = gauss_factor*gauss.pdf(d['distance'])

    if not np.any(mask):
        return False
    single_data = weights[mask]
    hist, bins = np.histogram(single_data['kappa'], bins=1000, range = (-0.2, 0.4))
    return (hist*gauss_factor/len(single_data), bins)



def get_bin_combos(bins, *args, **kwargs):
    bin_keys = [list(d.keys()) for d in bins.values()]
    combos = list(itertools.product(*bin_keys))
    return combos


def get_bins(dists: dict, data: dict, *args, **kwargs):
    bins = {}
    for name, dist_ in dists.items():
        fields = select_fields_by_weights(name, dist_, data, *args, **kwargs)
        bins.update({name: fields})

    return bins

def select_fields_by_weights(name, dist_, data, *args, **kwargs):
    nsigma = kwargs.get("nsigma", 2)
    binsize = kwargs.get("bin_size", 2)
    l_width = dist_[0] - dist_[1]
    r_width = dist_[2] - dist_[1]
    center = dist_[1]
    bins = np.arange(nsigma*l_width, nsigma*r_width+1, binsize)
    weights = data['ms_weights']
    med = weights[name].median()
    return_bins = {}
    for index, b in enumerate(bins):
        try:
            rb = bins[index+1]
        except IndexError:
            continue
        mask = (weights[name] > center + b) & (weights[name] < center + rb)
        if b < 0:
            f_distance = 0.5*(b+rb)/l_width
        else:
            f_distance = 0.5*(b+rb)/r_width
        return_bins.update({b: {'mask': mask, 'distance': f_distance}})
    

    return return_bins
        

def normalize_all(dists: dict, data: dict) -> dict:
    weights = data['ms_weights']
    median_n_gal = weights.gal.median()
    new_dists = {}
    for name, dist_ in dists.items():
        if name == 'gamma':
            w = weights['gamma']
    
            med = w.median()
            w = w*median_n_gal/med
            new_dists.update({name: [round(d*median_n_gal/med) for d in dist_]})
            weights[name] = w
        else:
            w = weights[name]
            med = w.median()
            w = median_n_gal*w/med
            new_dists.update({name: [round(d*median_n_gal) for d in dist_]})
            weights[name] = w

    return new_dists

def setup_ms(logger: logging.Logger, *args, **kwargs) -> Tuple:
    """
    Load relevant data from the millenium simulation
    """

    config = {'plane': kwargs.get("plane", False)}
    config['threads'] = kwargs.get("threads", 1)
    data = {}
    if not config['plane']:
        logger.critical("No redshift plane argument found.")
        raise KeyError("plane")

    config.update({'weight_path': kwargs.get("weight_path", False)})
    
    if not config['weight_path']:
        logging.critical("No path to load weights from!")
        raise KeyError("weight_path")
    
    try:
        config['weight_path'] = pathlib.Path(config['weight_path'])
    except TypeError:
        logging.critical("Path to weights cannot be cast as pathlib.Path")
        raise TypeError

    load_ms_weights(config['weight_path'], data)
    load_kappa_values(config, data, *args, **kwargs)

    config['gamma'] = kwargs.get("gamma", False)
    if config['gamma']:
        load_gamma_values(config, data)


    return config, data

def load_kappa_values(config: dict, data: dict, *args, **kwargs) -> None:
    plane = config['plane']
    directory = pathlib.Path(DATA_MANAGER.get_file_location({'datatype': 'kappa_maps', 'slice': str(plane)}))
    #PW: This is just the basename in my copy, may be worthwhile to make this more flexible
    cat = data['ms_weights']
    if 'kappa' in cat.columns:
        print("Using cached kappa values")
        return
    basename="GGL_los_8_{}_N"
    basepattern = "{}_{}"
    kappa_values = {}
    all_files = list(f.name for f in directory.glob('*.kappa'))
    all_kappa = np.zeros(len(cat), np.float32)
    point_lists = []
    file_list = []
    field_masks = []
    for i in range(8):
        for j in range(8):
            key = basepattern.format(i,j)
            id = basename.format(key)
            mask = (cat['field'] == key)
            points = list(zip(cat[mask].ra, cat[mask].dec))
            field_mask = data['ms_weights']['field'] == key
            field_cat = data['ms_weights'][field_mask]
            points = list(zip(field_cat.ra, field_cat.dec))
            fname = all_files[np.where([f.startswith(id) for f in all_files])[0][0]]
            fpath = directory / fname

            field_masks.append(field_mask)
            file_list.append(fpath)
            point_lists.append(points)

    nthreads = config['threads']
    with mp.Pool(nthreads) as p:
        kappas_ = p.map(get_kappa_values, list(zip(file_list, point_lists)))

    for index, gamma_vals in enumerate(kappas_):
        all_kappa[field_masks[index]] = gamma_vals
    
    data['ms_weights']['kappa'] = all_kappa
    rewrite_catalogs(data, *args, **kwargs)

def get_kappa_values(args):
    file = args[0]
    points = args[1]
    storage = np.zeros(len(points), np.float32)
    with open (file, 'rb') as d:
        f_data = np.fromfile(d, np.float32)
        f_data = np.reshape(f_data, (4096,4096))
        for index, point in enumerate(points):
            i,j = millenium_simulation.get_index_from_position(point[0]*u.deg, point[1]*u.deg)
            storage[index] = f_data[i,j]
    return storage


def load_gamma_values(config: dict, data: dict, *args, **kwargs) -> None:
    if "gamma" in data['ms_weights'].columns:
        print("Using cached gamma values")
        return
    
    plane = config['plane']
    gamma_dir = pathlib.Path(DATA_MANAGER.get_file_location({'datatype': 'gamma_maps', 'slice': str(plane)}))
    basename="GGL_los_8_{}_N"
    basepattern = "{}_{}"

    #PW: This is just the basename in my copy, may be worthwhile to make this more flexible
    gammas = np.zeros(len(data['ms_weights']), np.float32)
    all_files = list(gamma_dir.glob('*.gamma*'))
    file_lists = []
    point_lists = []
    field_masks = []

    for i in range(8):
        for j in range(8):
            key = basepattern.format(i,j)
            id = basename.format(key)
            fnames = list(filter(lambda x: x.name.startswith(id), all_files))
            if len(fnames) == 0:
                continue
            if len(fnames) != 2:
                raise FileNotFoundError("Found not enough or too many gamma files!")
            fpaths = [gamma_dir / fnames[0], gamma_dir / fnames[1]]
            field_mask = data['ms_weights']['field'] == key
            field_cat = data['ms_weights'][field_mask]
            points = list(zip(field_cat.ra, field_cat.dec))

            field_masks.append(field_mask)
            file_lists.append(fpaths)
            point_lists.append(points)
    
    nthreads = config['threads']
    with mp.Pool(nthreads) as p:
        gammas_ = p.map(get_gamma_values, list(zip(file_lists, point_lists)))

    for index, gamma_vals in enumerate(gammas_):
        gammas[field_masks[index]] = gamma_vals
    
    data['ms_weights']['gamma'] = gammas
    rewrite_catalogs(data, *args, **kwargs)

def rewrite_catalogs(data: dict, *args, **kwargs) -> None:
    cat = data['ms_weights']    
    cols = list(cat.columns)
    cols.remove('field')
    for field in cat['field'].unique():
        subcat = cat[cat['field'] == field]
        path = data['ms_files'][field]
        subcat.to_csv(path, index=False, columns=cols)
  
def get_gamma_values(vals, *args, **kwargs) -> np.array:
    """
    Gets the gamma values for a particular set of points in a particular subfield.
    This is slow, so written as it's own function for easier threading.

    Vals should be a list with these two entries (in order):
        gamma_file: list of pathlib.Paths. Should have length 2
        points: list of points
    Returns:
        gamma_values: list of gamma values at the points of interest
    
    """
    if len(vals) != 2:
        logging.error("Unable to unpack gamma values: wrong arguments!")
        return
    files = vals[0]
    points = vals[1]
    nfiles = len(files)
    if nfiles != 2:
        logging.error(f"Unable to unpack gamma values: Expected two files but got {nfiles}")
        return
    gamma_storage = np.zeros(len(points), np.float32)
    g1 = np.fromfile(files[0], np.float32)
    g2 = np.fromfile(files[1], np.float32)
    g1 = np.reshape(g1, (4096, 4096))
    g2 = np.reshape(g2, (4096, 4096))
    g = np.sqrt(g1**2 + g2**2)
    for index, pos in enumerate(points):
        ix, iy = millenium_simulation.get_index_from_position(pos[0]*u.deg, pos[1]*u.deg)
        g_val = g[ix, iy]
        gamma_storage[index] = g_val
    return gamma_storage



def load_ms_weights(path: pathlib.Path, data: dict, format="csv") -> None:
        """
        Loads weights output by the millenium simulation into a dataframe.
        Here, we just load everything into a single dataframe.

        Params:

        folder: location of the MS weights.

        """
        files = path.glob("*.csv")
        files = [f for f in files if not f.name.startswith('.')]
        ms_files = {}
        dfs = []

        for f in files:
            field = re.search(r"[0-9]_[0-9]", f.name)
            field = field.group()
            if format == "csv":
                f_path = path / f
                df = pd.read_csv(f_path)
                df['field'] = field
                dfs.append(df)
                ms_files.update({field: f})
        
        weights = pd.concat(dfs, ignore_index=True)
        if 'Unnamed: 0' in weights.columns:
            weights.drop(columns=['Unnamed: 0'], inplace=True)
        data.update({'ms_weights': weights, 'ms_files': ms_files})
def get_median(bins, hist):
    s = np.sum(hist)
    rs = 0.0
    for ix, v in enumerate(hist):
        rs += v
        if rs >= s/2:
            return bins[ix]

if __name__ == "__main__":

    centers1 = {'gal': 1.012, 'gamma':0.017}
    widths = {'gal': 0.05, 'gamma': 0.002}

    bins, hist = kappa_pdf( {'gal': [1.082, 1.112, 1.142], 'oneoverr': [1.08, 1.12, 1.16]}, plane=36, gamma=True, threads=8, weight_path = '/Users/patrick/Documents/Documents/Work/Research/LensEnv/ms/weighting2')
    med = get_median(bins, hist)
    import matplotlib.pyplot as plt
    plt.plot(bins[:-1], hist)
    plt.axvline(med, c='red')
    print(med)
    plt.show()
