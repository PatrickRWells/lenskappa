from asyncio import queues
from multiprocessing.sharedctypes import Value
from numpy.core.numeric import full
import pandas as pd
import os
import re
import logging
import itertools
import numpy as np
from pyparsing import col
from lenskappa.datasets.surveys.ms.ms import millenium_simulation
import random
import astropy.units as u
from scipy import stats
import math
import multiprocessing as mp
from lenskappa.utils import SurveyDataManager
import pathlib

def compute_histogram_range(combs, bins_, bin_range, normalized_weights, obs_centers, obs_widths, kappas, queue, tnum, *args, **kwargs):
        keys = list(bins_.keys())
        first = False
        thread_combs = combs[bin_range[0]: bin_range[1]]
        num = len(thread_combs)
        if num == 0:
            logging.critical("This thread recieved zero combinations to check!")
            raise ValueError

        logging.info("Thread {} got {} combinations".format(tnum, num))
        for index, comb in enumerate(thread_combs):
            if index%(num/10) == 0:
                logging.info("Thread {} completed {} of {} histograms".format(tnum, index, num))
            distance = []
            masks = []
            for i, key in enumerate(comb):
                bin = bins_[keys[i]][key]
                masks.append(bin['mask'])
                distance.append(bin['distance'])

            master_mask = np.all(masks, axis=0)
            if np.any(master_mask):
                if not first:
                    hist,bins = compute_single_histogram(normalized_weights, obs_centers, obs_widths, master_mask, distances=distance, kappas=kappas, *args, **kwargs)
                    first = True
                else:
                    dhist,bins = compute_single_histogram(normalized_weights, obs_centers, obs_widths, master_mask, distances=distance, kappas=kappas, *args, **kwargs)
                    hist += dhist
            else:
                pass
        logging.info("Thread {} finished work...".format(tnum))
        queue.put((bins, hist))

def compute_single_histogram(ms_weights, centers, widths, mask,
                            distances, kappas, min_kappa = -0.2, max_kappa = 1.0, kappa_bins=2000, *args, **kwargs):
        """
        Computes the historgram for a given n-dimensional bin output by get_bin_combos

        """

        gauss = stats.norm(0,1)
        #We weight the values based on their distance from the center of the
        #distribution.

        data = ms_weights[mask]
        #Get the weights that are in this bin



        fields = data['field'].unique()
        kappas = get_kappa_values(kappas, fields, *args, **kwargs)
        #Get the kappa values for all the fields (large 4deg x 4deg regions) found
        #in this subsample
        kappa_data = np.zeros(len(data))
        gauss_factor = 1.0
        for distance in distances:
            #Compute a gaussian factor for each weight being considered
            #Depending on its fractional distance from the center
            #PW: this should probably be moved to the get_bin_combos function
            #for clarity.
            gauss_factor *= gauss.pdf(distance)

        for val in data.field.unique():
            running_index = 0

            field_data = data[data.field == val]
            kappa_data_temp = np.zeros(len(field_data))
            for index,row in field_data.iterrows():
                indices = millenium_simulation.get_index_from_position(row.ra*u.deg, row.dec*u.deg)
                kappa_data_temp[running_index] = kappas[row.field][indices[0],indices[1]]
                running_index += 1

        running_index = 0
        for index, row in data.iterrows(): #Iterate over the fields being considered
            indices = millenium_simulation.get_index_from_position(row.ra*u.deg, row.dec*u.deg)
            #Find the index of the associated kappa point
            kappa_data[running_index] = kappas[row.field][indices[0],indices[1]]
            #Get the kappa value at that point
            running_index += 1
            for weight_name, weight_center in centers.items():
                weight_val = row[weight_name]
        #Create a histogram of the retrieved kappa values
        hist = np.histogram(kappa_data, bins=kappa_bins, range=(min_kappa, max_kappa))

        #Histogram is weighted by the distance from the center of the distribution
        #As well as the number of fields found
        return hist[0]*gauss_factor/len(data), hist[1]

def get_kappa_values(kappas, fields, *args, **kwargs):

    """
    Gets kappa values for a given subset of the 64 millenium simulation fields
    fields: list of (x,y) tuples, where x,y = [0,7]

    """

    outputs = {}
    for field in fields:
        outputs.update({field: kappas[field]})

    return outputs





class Kappa:
    def __init__(self) -> None:
        self._datamanager = SurveyDataManager("ms")


    def load_ms_weights(self, folder, format="csv", *args, **kwargs):
        """
        Loads weights output by the millenium simulation into a dataframe.
        Here, we just load everything into a single dataframe.

        Params:

        folder: location of the MS weights.

        """

        files = [f for f in os.listdir(folder) if f.endswith(format) and not f.startswith('.')]
        self._rootdir = folder
        self._ms_files = {}

        dfs = []
        for f in files:
            field = re.search(r"[0-9]_[0-9]", f)
            field = field.group()
            if format == "csv":
                path = os.path.join(folder, f)
                df = pd.read_csv(path)
                df['field']=field
                dfs.append(df)
                self._ms_files.update({field: f})
        
        self._weights = pd.concat(dfs, ignore_index=True)
        if 'Unnamed: 0' in self._weights.columns:
            self._weights.drop(columns=['Unnamed: 0'], inplace=True)

    def select_fields_by_weights(self, normalized_ms_weights, weight='gal', obs_center = 1.0, obs_width = 0.05, cwidth=2, bin_size=1):
        """
        Selects fields from the millenium simulation data, that have weights in the inputted range.
        This method only selects for a SINGLE weight.

        Params:

        normalized_ms_weights: Dataframe containing the normalized weights from the millenium simulation.
        weight: The name of the weight being used
        obs_center: Median weight value for the observed field of interest
        obs_width: Width of the distribution from the observed data.
        cwidth: How wide of a distribution to search.
            2 search from obs_center - obs_width to obs_center + obs_width
        bin_size: How large to make the search bins.

        """
        try:
            scaled_weight_vals = normalized_ms_weights[weight]
        except:
            logging.error("Unable to find weight {} in input data.".format(weight))
            return

        median_n_gals = normalized_ms_weights.gal.median()
        #This is a methods thing. We compare weight vals by first normalizing by the median
        #Then multiplying by the median number of galaxies in a MS field.
        #A bin width of 1 corresponds to increasing/decrease this median_n_gals by 1
        full_width = obs_width*median_n_gals*cwidth
        bin_width = bin_size
        nbins = int(2*cwidth/bin_size)
        center = median_n_gals*obs_center
        bins = [(center - cwidth)+i*bin_size for i in range(nbins+1)]

        #bins = [(center - full_width/2 + i*bin_width) for i in range(int(full_width/bin_width) + 1)]
        vals = {}
        bin_counts = np.zeros(len(bins))
        for ix, i in enumerate(bins):
            mask = (scaled_weight_vals >= i) & (scaled_weight_vals < i+bin_width)
            bin_counts[ix] = len(normalized_ms_weights[mask])
            vals.update({i: {'mask': mask, 'distance': (i-center)/center}})
            bin_counts[ix] = len(normalized_ms_weights[mask])
        return vals

    def get_bins(self, normalized_ms_weights, obs_centers, obs_widths, cwidth = 2, bin_size = 1, *args, **kwargs):

        """
        Gets all the bins given a set of weights, the normalized ms weights, and
        the distributions from the observed field.
        Number of weights scales as a^n, where n is the number of weights.

        Params:

        normalized_ms_weights: Normalized weights from the millenium simulation
        obs_centers: centers of the observed weight probability distributions
        obs_widths: widths of the observed weight probability distributions
        cwidth: How wide of a distribution to search.
            2 search from obs_center - obs_width to obs_center + obs_width
        """

        bins = {}
        for name, center in obs_centers.items():
            fields = self.select_fields_by_weights(normalized_ms_weights, name, center, obs_widths[name], cwidth, bin_size)
            bins.update({name: fields})
        return bins

    def get_bin_combos(self, bins, cwidth=4, bin_size=1, *args, **kwargs):
        """
        Combines the individual bins into combined, n-dimensional bins where n is the
        number of weights being considered.


        bins: Output of the get_bins function
        cwidth: see get_bins()
        bin_size see get_bins()
        """


        bin_keys = []
        for key, data in bins.items():
            bin_keys.append(list(data.keys()))


        combs = list(itertools.product(*bin_keys))
        return combs

    def get_kappa_values(self, fields, *args, **kwargs):

        """
        Gets kappa values for a given subset of the 64 millenium simulation fields
        fields: list of (x,y) tuples, where x,y = [0,7]

        """

        try:
            kappas = self._kappa_values
        except:

            self.load_kappa_values(*args, **kwargs)
        outputs = {}
        for field in fields:
            outputs.update({field: self._kappa_values[field]})

        return outputs


    def load_kappa_values(self, plane=36, *args, **kwargs):
        """
        Loads kappa values from disk

        Params:
        directory: Location of the kappa files
        """
        directory = pathlib.Path(self._datamanager.get_file_location({'datatype': 'kappa_maps', 'slice': str(plane)}))
        basename="GGL_los_8_{}_N"
        #PW: This is just the basename in my copy, may be worthwhile to make this more flexible
        basepattern = "{}_{}"
        self._kappa_values = {}
        all_files = list(f.name for f in directory.glob('*.kappa'))
        for i in range(8):
            for j in range(8):
                key = basepattern.format(i,j)
                id = basename.format(key)
                fname = all_files[np.where([f.startswith(id) for f in all_files])[0][0]]
                fpath = directory / fname
                with open (fpath, 'rb') as d:
                    data = np.fromfile(d, np.float32)
                    data = np.reshape(data, (4096,4096))
                    self._kappa_values.update({key: data})
    
    def attach_gamma(self, plane, *args, **kwargs):
        print("attaching gamma")
        catalog = self._weights
        gamma_directory = pathlib.Path(self._datamanager.get_file_location({'datatype': 'gamma_maps', 'slice': str(plane)}))
        basename="GGL_los_8_{}_N"
        #PW: This is just the basename in my copy, may be worthwhile to make this more flexible
        basepattern = "{}_{}"
        self._kappa_values = {}
        gammas = np.zeros(len(catalog), np.float32)
        all_files = list(f.name for f in gamma_directory.glob('*.gamma*'))
        file_lists = []
        point_lists = []
        field_masks = []
        for i in range(8):
            for j in range(8):
                key = basepattern.format(i,j)
                id = basename.format(key)
                fnames = list(filter(lambda x: x.startswith(id), all_files))
                if len(fnames) == 0:
                    continue
                if len(fnames) != 2:
                    raise FileNotFoundError("Found not enough or too many gamma files!")
                fpaths = [gamma_directory / fnames[0], gamma_directory / fnames[1]]
                field_mask = catalog['field'] == key
                field_cat = catalog[field_mask]
                points = list(zip(field_cat.ra, field_cat.dec))

                field_masks.append(field_mask)
                file_lists.append(fpaths)
                point_lists.append(points)
        
        print("getting gamma values!")
        nthreads = kwargs.get("nthreads", 1)
        with mp.Pool(nthreads) as p:
            print("Getting gammas")
            gammas_ = p.map(get_gamma_values, list(zip(file_lists, point_lists)))

        for index, gamma_vals in enumerate(gammas_):
            gammas[field_masks[index]] = gamma_vals
        
        catalog['gamma'] = gammas
        self._rewrite_catalogs()
    
    def _rewrite_catalogs(self, *args, **kwargs):
        cols = list(self._weights.columns)
        cols.remove('field')

        for field in self._weights['field'].unique():
            cat = self._weights[self._weights['field'] == field]
            path = pathlib.Path(self._rootdir) / self._ms_files[field]
            cat.to_csv(path, index=False, columns=cols)

    def normalize_weights(self, names, *args, **kwargs):
        """
        Noralized weights based so that the median value of a given
        weight = median_n_gal
        """
        median_n_gal = self._weights.gal.median()
        output=self._weights.copy(deep=True)
        for name in names:
            output[name] = median_n_gal*output[name]/output[name].median()
        return output

    def compute_histogram(self, ms_weights, centers, widths, mask,
                            distances, min_kappa = -0.2, max_kappa = 1.0, kappa_bins=2000, *args, **kwargs):
        """
        Computes the historgram for a given n-dimensional bin output by get_bin_combos

        """

        gauss = stats.norm(0,1)
        #We weight the values based on their distance from the center of the
        #distribution.

        data = ms_weights[mask]
        #Get the weights that are in this bin



        fields = data['field'].unique()
        kappas = self.get_kappa_values(fields, *args, **kwargs)
        #Get the kappa values for all the fields (large 4deg x 4deg regions) found
        #in this subsample
        kappa_data = np.zeros(len(data))
        gauss_factor = 1.0
        for distance in distances:
            #Compute a gaussian factor for each weight being considered
            #Depending on its fractional distance from the center
            #PW: this should probably be moved to the get_bin_combos function
            #for clarity.
            gauss_factor *= gauss.pdf(distance)

        for val in data.field.unique():
            running_index = 0

            field_data = data[data.field == val]
            kappa_data_temp = np.zeros(len(field_data))
            for index,row in field_data.iterrows():
                indices = millenium_simulation.get_index_from_position(row.ra*u.deg, row.dec*u.deg)
                kappa_data_temp[running_index] = kappas[row.field][indices[0],indices[1]]
                running_index += 1

        running_index = 0
        for index, row in data.iterrows(): #Iterate over the fields being considered
            indices = millenium_simulation.get_index_from_position(row.ra*u.deg, row.dec*u.deg)
            #Find the index of the associated kappa point
            kappa_data[running_index] = kappas[row.field][indices[0],indices[1]]
            #Get the kappa value at that point
            running_index += 1
            for weight_name, weight_center in centers.items():
                weight_val = row[weight_name]
        #Create a histogram of the retrieved kappa values
        hist = np.histogram(kappa_data, bins=kappa_bins, range=(min_kappa, max_kappa))

        #Histogram is weighted by the distance from the center of the distribution
        #As well as the number of fields found
        return hist[0]*gauss_factor/len(data), hist[1]

    def compute_kappa_pdf(self, obs_centers, obs_widths, plane=36, cwidth=2, bin_size = 1, min_kappa = -0.2, max_kappa = 1.0, kappa_bins=2000, nthreads=2, *args, **kwargs):
        #Compute and combine histograms for each bin into a single PDF for kappa
        print("Normalizing")

        if 'gamma' in obs_centers.keys():
            if 'gamma' not in self._weights.columns:
                self._weights['gamma'] = 0.0
                print("Attaching external shear values. This may take some time")
                print("These values will be stored to avoid having to do this in the future")
                self.attach_gamma(plane, nthreads=nthreads)
            gc = obs_centers['gamma']
            gw = obs_widths['gamma']
            med = np.median(self._weights['gamma'])
            if (med < 0):
                logging.warning("Median value of external shear is negative!")
            obs_centers['gamma'] = gc/med
            obs_widths['gamma'] = gw/med

        normalized_weights = self.normalize_weights(list(obs_centers.keys()))

        self.load_kappa_values(plane)
        print("Finding bins")
        bins_ = self.get_bins(normalized_weights, obs_centers=obs_centers, obs_widths=obs_widths, cwidth=cwidth, bin_size=bin_size, *args, **kwargs)
        keys = list(bins_.keys())
        print("Finding bin combos")
        combs = self.get_bin_combos(bins_, cwidth, bin_size, *args, **kwargs)
        #We shuffle the bin combinations, so that each thread gets (roughly) the same amount of work    
        random.shuffle(combs)

        first = False

        worker_threads = nthreads -1
        num_combs = len(combs)
        #replace this with a multithreaded version
        nperthread = math.floor(num_combs/worker_threads)
        thread_bins = [i*nperthread for i in range(worker_threads)]
        thread_bins.append(num_combs)
        processes = []
        queue = mp.Queue()
        print("Delegating to {} threads".format(worker_threads))
        for i, b in enumerate(thread_bins[:-1]):
            p = mp.Process(target=compute_histogram_range, args = (combs,bins_, (b, thread_bins[i+1]), normalized_weights, obs_centers, obs_widths, self._kappa_values,queue, i))
            p.start()
            processes.append(p)
            
        print("Getting histograms")
        hists = [queue.get() for _ in range(len(processes))]
        bins = hists[0][0]
        hist = hists[0][1]
        for item in hists[1:]:
            hist += item[1]

        for i,p in enumerate(processes):
            p.join()
            print("Thread {} joined!".format(i))



        print("Kappa pdf completed")
        return bins, hist           


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
        if index%10000 == 0:
            print(index)
        ix, iy = millenium_simulation.get_index_from_position(pos[0]*u.deg, pos[1]*u.deg)
        g_val = g[ix, iy]
        gamma_storage[index] = g_val
    return gamma_storage
