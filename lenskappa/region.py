from abc import ABC, abstractmethod
from astropy.coordinates.sky_coordinate import SkyCoord
from shapely import geometry
import logging
import astropy.units as u
import numpy as np
import time

from lenskappa.utils.multithreading import MultiThreadObject

class Region(ABC):
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
        own peril.
        Honestly, this should probably just be a subclass of a Shapely base object

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

        bounds = region.bounds
        ras = np.array([bounds[0], bounds[2]])
        decs = np.array([bounds[1], bounds[3]])
        if not ( (ras > 0.0).all() and (ras <= 360.0).all() ):
            logging.error("RAs must be between 0 and 360 degres, but got {}".format(ras))
            return None

        if not ( (decs >= -90.0).all() and (decs <= 90.0).all() ):
            logging.error("Declinations must be between -90 and 90 degres, but got {}".format(decs))
            return None



        super().__init__(center, region, *args, **kwargs)
        self._sample_type = "spherical"

    @property
    def skycoord(self):
        pass

    def is_ready(self, *args, **kwargs):
        try:
            sampler = self._sampler
        except:
            self._init_sampler(*args, **kwargs)
        return True

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
            logging.error("Attempted to draw a tile in the region, but the aperture does not have units")
            return False

        try:
            thread_num = kwargs['thread_num']
        except:
            thread_num = 0

        point_coords = self._draw_sample(thread_num)
        coord = SkyCoord(point_coords[0], point_coords[1], unit="deg")
        point = geometry.Point(point_coords)
        if not self._polygon.contains(point):
            #If the center of the tile falls outside the region, try again
            return self.generate_circular_tile(aperture=aperture, thread_num = thread_num, *args, **kwargs)

        return CircularSkyRegion(coord, aperture)


    def _init_sampler(self, *args, **kwargs):
        """
        Initialize the sampler for drawing regions.
        Currently, draws randomly off the surface of a sphere
        TODO: Test with iregularly-shaped region
        """
        if self._sample_type == 'spherical':
            self._init_spherical_sampler(*args, **kwargs)
        else:
            logging.error("Currently, only spherical sampling is implemented")

    def _init_spherical_sampler(self, *args, **kwargs):
        bounds = self._polygon.bounds
        ra1,ra2 = bounds[0], bounds[2]
        dec1, dec2 = bounds[1], bounds[3]
        ra_range = (min(ra1, ra2), max(ra1, ra2))
        dec_range = (min(dec1, dec2), max(dec1, dec2))
        phi_range = np.radians(ra_range)
        #Keeping everything in radians for simplicity
        #Convert from declination to standard spherical coordinates
        theta_range = (90. - dec_range[0], 90. - dec_range[1])
        #Area element on a sphere is dA = d(theta)d(cos[theta])
        #Sampling uniformly on the surface of a sphere means sampling uniformly
        #Over cos theta
        costheta_range = np.cos(np.radians(theta_range))

        self._low_sampler_range = [phi_range[0], costheta_range[0]]
        self._high_sampler_range = [phi_range[1], costheta_range[1]]
        self._sampler = MultiThreadObject(np.random.default_rng)

    def _draw_sample(self, thread_num = 0):
        try:
            sampler = self._sampler
        except:
            self._init_sampler()

        if self._sample_type == "spherical":
            vals = self._sampler[thread_num].uniform(self._low_sampler_range, self._high_sampler_range)
            ra = np.degrees(vals[0])
            theta = np.degrees(np.arccos(vals[1]))
            dec = 90 - theta
            return (ra, dec)

    def contains(self, obj):
        if type(obj) == CircularSkyRegion:
            return self._polygon.contains(obj._polygon)
        elif type(obj) == geometry.Point:
            return self._polygon.contains(obj)

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
