from abc import ABCMeta, abstractmethod
import logging
from astropy.coordinates.sky_coordinate import SkyCoord
import astropy.units as u
import numpy as np
import math
from numpy.lib.arraysetops import isin
from numpy.lib.polynomial import poly
from shapely import geometry
import time
import matplotlib.pyplot as plt
from shapely.geometry import geo

class Region(metaclass=ABCMeta):
    def __init__(self, center, polygon, *args, **kwargs):
        self._center = center
        self._polygon = polygon

    def get_polygon(self):
        return self._polygon

    def add_subregion(self, name, center, input_region, override = False, *args, **kwargs):
        """
        Adds a subreegion to the current region.
        This class also supports adding sub-subregions through the "propogate_subregion" method
        Subregions must be fuly enclosed by their parent regions. This can be overriden, but do so at your
        own peril
        Parameters:
            name: Name for the subregion. Can be anything that works as a dictionary key
            center: Tuple with (x,y) coordinates of the geometric center or the region
            polygon: Polygon or list of (x,y) points defining the corners
            override: Whether to override the interior requirement (see ablove)
        Returns:
            True if the subregion was sucessfully added, false otherwise
        """
        try:
            subregions = self._subregions
        except:
            self._subregion_centers = {}
            self._subregions = {}

        within = self._polygon.contains(input_region)
        if not within and not override:
            logging.error("Tried to initialize subregion {}, but it is not contained with its parent region".format(name))
            return False
        elif isinstance(input_region, Region):
            self._subregions.update({name: input_region})
        else:
            self._subregions.update({name: self.build_region(center, input_region, *args, **kwargs)})
        self._subregion_centers.update({name: center})
        return True
    
    @property
    def area(self):
        return self._polygon.area
    
    @property
    def center(self):
        return self._center
    
    def get_type(self):
        poly_type = type(self._polygon)
        if poly_type == geometry.Polygon:
            return "polygon"
        elif poly_type == geometry.Point:
            return "circle"
        else:
            logging.error("Invalid polygon type encountered: {}".format(poly_type))


    @staticmethod
    def build_region(center, input_region, *args, **kwargs):
        pass

    @abstractmethod
    def intersects(self, second_region, *args, **kwargs):
        pass

    def intersects(self, second_region, *args, **kwargs):
        try:
            reg = second_region.get_polygon()
        except:
            reg = second_region
        return self._polygon.intersects(reg)
    
    def get_subregion_intersections(self, second_region):
        sub = []
        for name, subregion in self._subregions.items():
            if subregion.intersects(second_region):
                sub.append(name)
        return sub
    
    @classmethod
    def union(cls, regions, *args, **kwargs):
        if type(regions) != list:
            logging.error("Unable to combine regions, expected a list")
            return None
        shapely_objs = [region._polygon for region in regions]
        from shapely.ops import unary_union
        combined = unary_union(shapely_objs)
        return cls(combined.centroid, combined)

        


class SkyRegion(Region):
    def __init__(self, center, region, *args, **kwargs):
        super().__init__(center, region, *args, **kwargs)
    
    @property
    def skycoord(self):
        pass

    def build_region(self, center, input_region, *args, **kwargs):
        return SkyRegion(center, input_region)

    def generate_circular_tile(self, aperture, *args, **kwargs):
        """
        Generates a circular tile, somewhere in the region.
        In the past this was done by fully tiling the region before doing the weights
        Doing this randomly is advantageous for several reasons:
            1. Handles non-rectangular regions natively (just rejects samples that fall outside)
            2. Makes the code easier to parallelize
            3. Allows for easier testing
        Parameters:
            aperture: <astropy.Quantity> radius of the circle
        returns:
            region: Region objects defining the tile


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
        coord = SkyCoord(point_coords[0], point_coords[1], unit="deg")
        point = geometry.Point(point_coords)
        if not self._polygon.contains(point):
            #If the center of the tile falls outside the region, try again
            return self.generate_circular_tile(aperture=aperture)
        
        return CircularSkyRegion(coord, aperture)
        
        
    def _init_sampler(self):
        """
        Initialize the sampler for drawing regions
        TODO: Test with iregularly-shaped region
        """
        bounds = self._polygon.bounds
        x1, x2 = bounds[0], bounds[2]
        x_range = (min(x1, x2), max(x1, x2))
        y1, y2 = bounds[1], bounds[3]
        y_range = (min(y1, y2), max(y1, y2))

        self._low_sampler_range = [x_range[0], y_range[0]]
        self._high_sampler_range = [x_range[1], y_range[1]]

        self._sampler = np.random.default_rng(int(time.time()))

    
    def contains(self, point):
        return self._polygon.contains(point)
    
class CircularSkyRegion(SkyRegion):

    def __init__(self, center, radius, *args, **kwargs):
        """
        Named "CircularSkyRegion to avoid conflict with the astropy.Regions object
        "CircleSkyRegion
        """
        try:
            self._coord = center
            x = center.ra.degree
            y = center.dec.degree
            point = geometry.Point(x,y)
        except:
            logging.error("Unable to initialize point {}. Expected a SkyCoord".format(center))
            return(None)
        
        try:
            rad = radius.to(u.degree)
            self._radius = rad
            circle = point.buffer(rad.value)
        except:
            logging.error("Unable creat circle of radius {}. Radius should be an astropy quantity".format(radius))
        super().__init__(center, circle)
    
    @property
    def skycoord(self):
        return self._coord, self._radius
