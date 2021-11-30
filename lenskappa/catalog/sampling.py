from abc import ABC, abstractmethod
import numpy as np
import logging
from copy import deepcopy

from lenskappa.catalog import params

class Distribution(ABC):


    def __init__(self, dist_type, *args, **kwargs):
        self._type = dist_type
    
    @abstractmethod
    def generate_samples(self, n, *args, **kwargs):
        pass
    
    def get_sample(self, *args, **kwargs):
        for sample in np.transpose(self._samples):
            yield sample
    
    def get_all_samples(self, *args, **kwargs):
        try:
            return self._samples
        except:
            logging.error("No samples have been generated!")
            return

    def apply_boolean_mask(self, mask):
        cp = deepcopy(self)
        cp._samples = cp._samples[mask]
        return cp

class ArbitraryDistribution(Distribution):
    def __init__(self, samples, *args, **kwargs):
        super().__init__("arbitrary")
        self._samples = samples
    def generate_samples(self, n, *args, **kwargs):
        pass

class GaussianDistribution(Distribution):
    def __init__(self, centers, widths, *args, **kwargs):
        super().__init__("gaussian", *args, **kwargs)
        self._centers = centers
        self._widths = widths

    @property
    def num_samples(self):
        return self._samples.shape[1]

    def generate_samples(self, n, *args, **kwargs):
        num_objs = len(self._centers)
        rng = np.random.default_rng()
        samples = np.zeros(shape=(num_objs,n), dtype=np.float64)
        scaled_samples = np.zeros(shape=(num_objs,n), dtype=np.float64)
        rng.standard_normal(size=(num_objs,n), out=samples)
        for index, center in enumerate(self._centers):
            scaled_samples[index] = self._widths[index]*samples[index] + center

        self._samples = scaled_samples
    
def process_samples_from_array(catalog, id_column, samples, parname):
    n_samples = len(samples.columns)-2 #Assume first column contains object IDs
    data = np.zeros((len(catalog), n_samples))
    ids = samples['objid']
    for index, row in samples.iterrows():
        sample_id = ids[index]
        mask =  np.where(catalog[id_column] == sample_id)
        row_i = mask[0][0]
        data[row_i] = row.drop(labels=['Unnamed: 0', 'objid'])
    dist = ArbitraryDistribution(data)
    catalog.attach_samples(parname, dist)