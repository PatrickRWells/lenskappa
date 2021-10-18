import pandas as pd
import os
import re
import logging
import itertools
import numpy as np
from lenskappa.surveys.ms.ms import millenium_simulation
import astropy.units as u
from scipy import stats

class Kappa:
    def __init__(self) -> None:
        pass

    def load_ms_weights(self, folder, format="csv", *args, **kwargs):
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

    def select_fields_by_weights(self, normalized_ms_weights, weight='gal', obs_center = 1.0, obs_width = 0.05, cwidth=2, bin_size=1):
        try:
            scaled_weight_vals = normalized_ms_weights[weight]
        except:
            logging.error("Unable to find weight {} in input data.".format(weight))
            return 
        #Normalize weight vals
        median_n_gals = normalized_ms_weights.gal.median()

        full_width = obs_width*median_n_gals*cwidth
        bin_width = bin_size

        center = median_n_gals*obs_center
        import matplotlib.pyplot as plt 
        bins = [(center - full_width/2 + i*bin_width) for i in range(int(full_width/bin_width) + 1)]
        vals = {}
        for i in bins:
            mask = (scaled_weight_vals >= i) & (scaled_weight_vals < i+bin_width)
            vals.update({i: {'mask': mask, 'distance': (i-center)/center}})
        return vals

    def get_bins(self, normalized_ms_weights, obs_centers, obs_widths, cwidth = 2, bin_size = 1, *args, **kwargs):

        bins = {}
        for name, center in obs_centers.items():
            fields = self.select_fields_by_weights(normalized_ms_weights, name, center, obs_widths[name], cwidth, bin_size)
            bins.update({name: fields})
        return bins

    def get_bin_combos(self, bins, cwidth=4, bin_size=1, *args, **kwargs):
        bin_keys = []
        for key, data in bins.items():
            bin_keys.append(list(data.keys()))


        combs = list(itertools.product(*bin_keys))
        keys = list(bins.keys())
        output = {}
        print("Getting bin combinations!")
        for index, c in enumerate(combs):
            distances = np.zeros(len(keys))
            master_mask = np.ones(len(self._weights), dtype=bool)
            for index, val in enumerate(c):
                mask = bins[keys[index]][val]['mask']
                master_mask = master_mask&mask
                distances[index] = bins[keys[index]][val]['distance']
            output.update({c: {'mask': master_mask, 'distances': distances}})
        return output

    def get_kappa_values(self, fields, *args, **kwargs):
        try:
            kappas = self._kappa_values
        except:
            
            self.load_kappa_values(*args, **kwargs)
        
        outputs = {}
        for field in fields:
            outputs.update({field: self._kappa_values[field]})

        return outputs


    def load_kappa_values(self, directory, *args, **kwargs):
        basename="GGL_los_8_{}_N"
        basepattern = "{}_{}"
        self._kappa_values = {}
        all_files = [f for f in os.listdir(directory) if not f.startswith('.')]
        for i in range(8):
            for j in range(8):
                key = basepattern.format(i,j)
                id = basename.format(key)
                fname = all_files[np.where(id)[0][0]]
                fpath = os.path.join(directory, fname)
                with open (fpath, 'rb') as d:
                    data = np.fromfile(d, np.float32)
                    data = np.reshape(data, (4096,4096))
                    self._kappa_values.update({key: data})
             
    def normalize_weights(self, names, *args, **kwargs):
        median_n_gal = self._weights.gal.median()
        output=self._weights.copy(deep=True)
        for name in names:
            output[name] = median_n_gal*output[name]/output[name].median()
        return output

    def compute_histogram(self, ms_weights, centers, widths, mask, distances, min_kappa = -0.2, max_kappa = 1.0, kappa_bins=2000, *args, **kwargs):
        gauss = stats.norm(0,1)
        data = ms_weights[mask]
        fields = data[mask]['field'].unique()
        kappas = self.get_kappa_values(fields, *args, **kwargs)
        kappa_data = np.zeros(len(data))
        gauss_factor = 1.0
        for distance in distances:
            gauss_factor *= gauss.pdf(distance)
        running_index = 0
        for index, row in data.iterrows():
            indices = millenium_simulation.get_index_from_position(row.ra*u.deg, row.dec*u.deg)
            kappa_data[running_index] = kappas[row.field][indices[0],indices[1]]
            running_index += 1
            for weight_name, weight_center in centers.items():
                weight_val = row[weight_name]

        hist = np.histogram(kappa_data, bins=kappa_bins, range=(min_kappa, max_kappa))
        return hist[0]*gauss_factor/len(data), hist[1]

    def compute_histograms(self, obs_centers, obs_widths, cwidth=2, bin_size = 1, min_kappa = -0.2, max_kappa = 1.0, kappa_bins=2000, *args, **kwargs):
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






