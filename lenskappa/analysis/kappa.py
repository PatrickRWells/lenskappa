from numpy.core.numeric import full
import pandas as pd
import os
import re
import logging
import itertools
import numpy as np
from lenskappa.datasets.surveys.ms.ms import millenium_simulation
import astropy.units as u
from scipy import stats

class Kappa:
    def __init__(self) -> None:
        pass

    def load_ms_weights(self, folder, format="csv", *args, **kwargs):
        """
        Loads weights output by the millenium simulation into a dataframe.
        Here, we just load everything into a single dataframe.

        Params:

        folder: location of the MS weights.

        """

        files = [f for f in os.listdir(folder) if f.endswith(format)]
        dfs = []
        for f in files:
            field = re.search(r"[0-9]_[0-9]", f)
            field = field.group()
            if format == "csv":
                path = os.path.join(folder, f)
                df = pd.read_csv(path)
                df['field']=field
                dfs.append(df)
        self._weights = pd.concat(dfs, ignore_index=True)
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

        center = median_n_gals*obs_center

        bins = [(center - full_width/2 + i*bin_width) for i in range(int(full_width/bin_width) + 1)]

        vals = {}
        bin_counts = np.zeros(len(bins))
        for ix, i in enumerate(bins):
            mask = (scaled_weight_vals >= i) & (scaled_weight_vals < i+bin_width)
            bin_counts[ix] = len(normalized_ms_weights[mask])
            vals.update({i: {'mask': mask, 'distance': (i-center)/center}})
            bin_counts[ix] = len(normalized_ms_weights[mask])
            import matplotlib.pyplot as plt
        return vals

    def get_bins(self, normalized_ms_weights, obs_centers, obs_widths, cwidth = 2, bin_size = 1, *args, **kwargs):

        """
        Gets all the bins given an set of weights, the normalized ms weights, and 
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
        keys = list(bins.keys())
        output = {}
        combo_counts = np.zeros(len(combs))
        for index_1, c in enumerate(combs):
            distances = np.zeros(len(keys))
            master_mask = np.ones(len(self._weights), dtype=bool)
            for index_2, val in enumerate(c):
                mask = bins[keys[index_2]][val]['mask']
                master_mask = master_mask&mask
                distances[index_2] = bins[keys[index_2]][val]['distance']
            output.update({c: {'mask': master_mask, 'distances': distances}})
        return output

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


    def load_kappa_values(self, directory, *args, **kwargs):
        """
        Loads kappa values from disk

        Params:
        directory: Location of the kappa files
        """
        basename="GGL_los_8_{}_N"
        #PW: This is just the basename in my copy, may be worthwhile to make this more flexible
        basepattern = "{}_{}"
        self._kappa_values = {}
        all_files = [f for f in os.listdir(directory) if not f.startswith('.')]
        for i in range(8):
            for j in range(8):
                key = basepattern.format(i,j)
                id = basename.format(key)
                fname = all_files[np.where([f.startswith(id) for f in all_files])[0][0]]
                fpath = os.path.join(directory, fname)
                with open (fpath, 'rb') as d:
                    data = np.fromfile(d, np.float32)
                    data = np.reshape(data, (4096,4096))
                    self._kappa_values.update({key: data})
             
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
        


        fields = data[mask]['field'].unique()
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
        import matplotlib.pyplot as plt

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

    def compute_kappa_pdf(self, obs_centers, obs_widths, cwidth=2, bin_size = 1, min_kappa = -0.2, max_kappa = 1.0, kappa_bins=2000, *args, **kwargs):
        #Compute and combine histograms for each bin into a single PDF for kappa
        normalized_weights = self.normalize_weights(list(obs_centers.keys()))
        bins = self.get_bins(normalized_weights, obs_centers=obs_centers, obs_widths=obs_widths, cwidth=cwidth, bin_size=bin_size, *args, **kwargs)
        combs = self.get_bin_combos(bins, cwidth, bin_size, *args, **kwargs)
        print("Number of combinations to check {}".format(len(combs)))
        master_mask = np.zeros(len(self._weights), bool)
        first = False
        num_combs = len(combs)
        for index, (bin, data) in enumerate(combs.items()):
            if index%100 == 0:
                print("Completed {} out of {} histograms".format(index, num_combs))
            mask = data['mask']
            if np.any(mask):
                if not first:
                    hist,bins = self.compute_histogram(normalized_weights, obs_centers, obs_widths, mask, distances=data['distances'], *args, **kwargs)
                    first = True
                else:
                    dhist,bins = self.compute_histogram(normalized_weights, obs_centers, obs_widths, mask, distances=data['distances'], *args, **kwargs)
                    hist += dhist
            else:
                pass

        return bins, hist





