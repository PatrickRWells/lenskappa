from abc import ABC, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
import logging

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


class GaussianDistribution(Distribution):
    def __init__(self, centers, widths, *args, **kwargs):
        super().__init__("gaussian", *args, **kwargs)
        self._centers = centers
        self._widths = widths

    def apply_boolean_mask(self, mask):
        output =  GaussianDistribution(self._centers, self._widths)
        output._samples = self._samples[mask]
        return output

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
    


