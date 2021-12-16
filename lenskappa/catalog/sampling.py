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

class DistributionArray(ABC):

    def __init__(self, *args, **kwargs):
        """
        A distribution array is used to attach a collection of distributions
        to a a set of particular values of parameters in a catalog.
        Originally, this was implemented to allow us to create artificial
        pdfs for redshifts for an object in a catalog during weighting.
        """
    pass


class GaussianDistributionArray(DistributionArray):
    def __init__(self, *args, **kwargs):
        """
        A distribution array, where each distribution is a gaussian
        """
        pass

    def generate_grid(self, params, limits, bins, *args, **kwargs):
        """
        Generate a n-dimensional grid of distributions, where n
        is the number of parameters.
        """
        length = len(params)
        if not all(len(l) == length for l in [limits, bins]):
            logging.error("Could not init distributions object"\
                          "Expected params, size, and limits to be the same length")
            return
        
        if not all(len(l) == 2 for l in limits):
            logging.error("Unable to init distributions object"\
                "Expected two limits for each parameter")
            return

        self._axes = {p: np.linspace(*limits[i], bins[i]) for i,p in enumerate(params)}

    def _declare_dists(self, *args, **kwargs):
        shape = [len(coord) for coord in self._axes.values()]
        self._distributions = np.empty(shape=shape, dtype=tuple)
        return self._distributions

    def add_distribution(self, coordinates, center, width):
        """
        Attaches a distribution to a given location in parameter space.
        Center and width are the center and width of the gaussian
        """
        try:
            dists = self._distributions
        except:
            dists = self._declare_dists()

        indices = np.zeros(len(self._axes), dtype=int)
        for i, (coordinate, axis) in enumerate(self._axes.items()):
            try:
                coord_val = coordinates[coordinate]
            except:
                logging.error("Unable to attach sample")
                return
            index = np.searchsorted(axis, coord_val, "left")
            indices[i] = index
        dists[tuple(indices)] = (center, width)
    
    def get_distribution(self, coordinates):

        """
        Return the distribution closest to the parameter values passed
        Coordinates: {parameter_name: parameter_value}
        """
        indices = np.zeros(len(self._axes), dtype=int)
        for i, (coordinate, axis) in enumerate(self._axes.items()):
            try:
                coord_val = coordinates[coordinate]
            except:
                logging.error("Unable to attach sample")
                return
            index = np.searchsorted(axis, coord_val, "left")
            indices[i] = index
        
        return self._distributions[tuple(indices)]


    def draw(self, coordinates, n):
        """
        Return n samples from the distribution closest to the parameter values passed
        Coordinates: {parameter_name: parameter_value}
        n: number of samples to draw. 
        """

        try:
            rng = self._rng
        except:
            self._rng = np.random.default_rng()
        
        dist = self.get_distribution(coordinates)
        if dist is not None:
            vals = self._rng.standard_normal(size=n)
            samples = dist[0] + vals*dist[1]
            return samples
        else:
            logging.error("Distribution at point {} not initialized".format(coordinates))

    def draw_many(self, coordinates, n):
        """
        Draw samples from many different distributions
        coordinates: {parameter_name: [parameter_values]...}
        n: number of samples to draw from each distribution

        Coordinate arrays should be the same length
            coordinates will be paired element-wise
        """
        coordinate_values = [coordinates[c] for c in self._axes]
        coordinate_len = len(coordinate_values[0])
        if not all(len(cv) == coordinate_len for cv in coordinate_values):
            logging.error("When drawing from many distributions, all"\
                            " coordinate arrays must be the same length") 
        k = list(self._axes.keys())

        value_dicts = list(map(lambda x: {a: x[i] for i, a in enumerate(k)}, (zip(*coordinate_values))))
        #What a terrible line amiright?
        distributions = [self.get_distribution(d) for d in value_dicts]

        try:
            rng = self._rng
        except:
            self._rng = np.random.default_rng()

        normalized_values = self._rng.standard_normal((coordinate_len, n))
        values = np.zeros((coordinate_len, n), dtype=float)
        for index, row in enumerate(normalized_values):
            dist = distributions[index]
            values[index] = dist[0] + row*dist[1]

        return values


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
