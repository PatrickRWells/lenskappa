from abc import ABCMeta, abstractmethod
import logging
import astropy.units as u
import numpy as np
import math
from numpy.lib.polynomial import poly
from shapely import geometry
import time
import matplotlib.pyplot as plt
from shapely.geometry import geo

class Region(metaclass=ABCMeta):
    def __init__(self, center, region, *args, **kwargs):
        self._center = center
        self._region = region
    
    def add_subregion(self, name, center, input_region, override = False, *args, **kwargs):
        """
        Adds a subreegion to the current region.
        This class also supports adding sub-subregions through the "propogate_subregion" method
        Subregions must be fuly enclosed by their parent regions. This can be overriden, but do so at your
        own peril
        Parameters:
            name: Name for the subregion. Can be anything that works as a dictionary key
            center: Tuple with (x,y) coordinates of the center
            input_region: Polygon or list of (x,y) points defining the corners
            override: Whether to override the interior requirement (see ablove)
        Returns:
            True if the subregion was sucessfully added, false otherwise
        """
        try:
            subregions = self._subregions
        except:
            self._subregion_centers = {}
            self._subregions = {}

        within = self._region.contains(input_region)
        if not within and not override:
            logging.error("Tried to initialize subregion {}, but it is not contained with its parent region".format(name))
            return False
        elif isinstance(input_region, Region):
            self._subregions.update({name: input_region})
        else:
            self._subregions.update({name: self.build_region(center, input_region, *args, **kwargs)})
        self._subregion_centers.update({name: center})
        return True
    
            
    @staticmethod
    def build_region(center, input_region, *args, **kwargs):
        pass

    @abstractmethod
    def overlaps(self, second_region, *args, **kwargs):
        pass

class SkyRegion(Region):
    def __init__(self, center, region, *args, **kwargs):
        super().__init__(center, region, *args, **kwargs)

    def build_region(self, center, input_region, *args, **kwargs):
        return SkyRegion(center, input_region)

    def overlaps(self, second_region):
        return False

    def generate_tile(self, aperture, *args, **kwargs):
        """
        Generates a circular tile, somewhere in the region.
        In the past this was done by fully tiling the region before doing the weights
        Doing this randomly is advantageous for several reasons:
            1. Handles non-rectangular regions natively (just rejects sampels that fall outside)
            2. Makes the code easier to parallelize
            3. Allows for easier testing
        Parameters:
            aperture: <astropy.Quantity> radius of the circle
        returns:
            region: Region objects defining the tile

        
        TODO: Implement different tiling methods
        TODO: Implement for non-rectangular region

        """
        try:
            size = aperture.to(u.deg)
        except:
            logging.error("Attempted to tile the region, but the aperture does not have units")
            return False
        
        try:
            sampler = self._sampler
        except:
            self._init_sampler()
            sampler = self._sampler
        
        point_coords = sampler.uniform(self._low_sampler_range, self._high_sampler_range)
        point = geometry.Point(point_coords)
        if not self._region.contains(point):
            #If the center of the tile falls outside the region, try again
            return self.generate_tile(aperture=aperture)
        
        return point.buffer(size.value)
        
        
    def _init_sampler(self):
        """
        Initialize the sampler for drawing regions
        TODO: Test with iregularly-shaped region
        """
        bounds = self._region.bounds
        x1, x2 = bounds[0], bounds[2]
        x_range = (min(x1, x2), max(x1, x2))
        y1, y2 = bounds[1], bounds[3]
        y_range = (min(y1, y2), max(y1, y2))

        self._low_sampler_range = [x_range[0], y_range[0]]
        self._high_sampler_range = [x_range[1], y_range[1]]

        self._sampler = np.random.default_rng(int(time.time()))