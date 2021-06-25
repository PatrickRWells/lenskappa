from astropy.coordinates.sky_coordinate import SkyCoord
import pandas as pd
import numpy as np
from astropy import wcs, coordinates, units
import astropy.units as u
from shapely import geometry
from copy import copy
import logging
import shapely
from shapely.geometry import geo
import toml
from abc import ABCMeta, abstractmethod

def require_points(fn):
    """
    Decorator to check that the object points in a catalog have
    been initialized
    """
    def wrapper(self, *args, **kwargs):
        if not self._points_initialized:
            self._init_points()
        return fn(self, *args, **kwargs)
    return wrapper


def require_validation(params, *args, **kwargs):
    """
    Decorator to ensure a particular set of parameters
    Have been validate against the input catalog
    eg:
    @require_validation(['x', 'y'])
    def func(self, *args, **kwargs)

    This would ensure that the parameters x, and y have been mapped
    to the appropriate catalog columns and their types are valid

    """
    def outer(fn, *args, **kwargs):

        def wrapper(self, *args, **kwargs):
            try:
                checkvalid = self._valid
            except:
                logging.error("Catalog parameters have not been validated!")
                return None
            isvalid = [checkvalid[name] for name in params]
            if np.all(isvalid):
                return fn(self, *args, **kwargs)
            else:
                logging.error("Catalog parameter(s) {} are not valid".format([par for index, par in enumerate(pars) if not isvalid[index]]))
        return wrapper
    return outer

class Catalog(metaclass=ABCMeta):

    def __init__(self, cat, parmap = {}, partypes = {}, *args, **kwargs):
        self._cat = cat
        self._needed_parameter_types = partypes
        self._parmap = parmap
        if parmap:
            self._validate_parmap()

    def __getitem__(self, key):
        """
        Allows for masking as with a usual dataframe
        """
        return self._cat[key]

    def __setitem__(self, key, data):
        """
        Allows for setting values as with a usual dataframe
        """
        self._cat[key] = data

    def __len__(self):
        """
        Returns the length of the underlying data
        """
        return len(self._cat)

    def __str__(self):
        return self._cat.__str__()


    @classmethod
    def read_csv(cls, file, config = None, *args, **kwargs):
        """Read in a file from a CSV
        Arguments:
            - file (str): path to CSV file
            - config (dict): dictionary of config parameters
        Rerturns:
            catalog.Catalog
        """
        logging.info("Reading in file {}".format(file))
        try:
            cat = pd.read_csv(file)
            return cls(cat, *args, **kwargs)
        except:
            logging.error("Pandas could not read file {}".format(file))
            return None

    def load_params(self, input_map, *args, **kwargs):
        """
        Used for marking catalog columns with standard labels
        e.g. The code expects columns "ra" and "dec", while
        the catalog is labeled with "RaMean" and "DecMean"

        Parameters:
            input_map {str: str} - key-value pairs where the key
                is the standard column name and the value is the
                actual column name.
        """

        try:
            parmap = self._parmap
        except:
            self._parmap = {}

        for par, parmap in self._parmap.items():
            try:
                input_parmap = input_map[par]
                self._parmap.update({par: input_parmap})
            except:
                logging.warning("No rename found for parameter {}, defaulting to {}".format(par, parmap))
        self._validate_parmap()

    def _validate_parmap(self):
        """
        Checks to make sure that all necessary parameters are present
        In the catalog and the values are of the expected type.
        Note: Currently checks for type by attempting to cast the column as
        the expected type. This should handle problems like getting a
        np.float instead of float, but may have unintended side effects
        #TODO Implement re-evaluation when parameters are changed
        """
        self._valid = {}
        for par, partype in self._needed_parameter_types.items():

            single_parmap = self._parmap[par]
            try:
                col = self._cat[single_parmap]
            except:
                logging.warning("Unable to find column {} in data".format(single_parmap))
                self._valid.update({par: False})
                continue

            try:
                coldata = col.astype(partype)
                self._valid.update({par : True})
            except:
                logging.error("Unable to cast column {} as type {}".format(single_parmap, partype))
                self._valid.update({par : False})


    @abstractmethod
    def add_subregion(self, *args, **kwargs):
        pass

    @abstractmethod
    def filter_by_subregion(self, *args, **kwargs):
        pass
    
    @abstractmethod
    def get_objects_in_region(self, region):
        pass

class Catalog2D(Catalog):
    def __init__(self, cat, *args, **kwargs):
        needed_parameter_types = {'x': float, 'y': float}
        super().__init__(cat, partypes=needed_parameter_types, *args, **kwargs)
        self._points_initialized = False
        self._subregions = {}

    def add_subregion(self, name, polygon, *args, **kwargs):
        """
        Adds a subregion to the catalog
        This is mostly just a convinient way of subdividing catalog data

        params:
            name <str>: Name for the subregion
            polygon <shapely.geometry.Polygon>: Polygon defining the region

        """
        polytype = type(polygon)
        if polytype in [geometry.Polygon, geometry.Point]:
            self._subregions.update({name: {'poly': polygon, 'type': polytype } } )
        else:
            logging.error("Tried to add a subregion, but the region was not a polygon")

    @require_validation(params=['x', 'y'])
    def filter_by_subregion(self, name, *args, **kwargs):
        """
        Filters data by a particular subregion.
        This method is lazy, in that it doesn't actually
        evaluate which objects are in the subregion until this method
        is called for the first time

        """
        if name not in self._subregions.keys():
            logging.error("Subregion {} does not exist".format(name))
            return

        if not self._subregions[name]['mask']:
            self._init_subregion(name)
        return self._subregions[name]['mask']

    @require_validation(params=['x', 'y'])
    def _init_points(self):
        """
        Initialize a set of shapely points and store them in the catalog
        TODO: Stress-test with millions of points
        """
        x = self._cat[self._parmap['x']]
        y = self._cat[self._parmap['y']]
        self._points = [geometry.Point(xy) for xy in list(zip(x,y))]
        self._cat['mask'] = self._points
        self._points_initialized = True

    @require_validation(params=['x', 'y'])
    def _init_subregion(self, name):
        """
        Initializes subregion by finding which points fall inside of it
        Results are stored as a mask that can be applied to _cat
        Method assumes that _init_points has already been run
        """
        poly = self._subregions[name]['poly']
        #Type validation of polygons is handled when the subregion is created
        mask = [point.within(poly) for point in self._points]
        initial_len = len(self._cat)
        final_len = len(self._cat[mask])
        logging.info("Region {} masks out {} objects".format(name, initial_len-final_len))
        self._subregions[name]['mask'] = mask
    
    def get_objects_in_region(self, region):
        #Note, this is a very simple implementationt that will be inefficient for very large regions
        #It is recommended that you override this for larger catalogs 
        mask = np.array([region.contains(point) for point in self._points])
        return self[mask]
        


class SkyCatalog2D(Catalog2D):

    def __init__(self, cat, *args, **kwargs):
        parmap = ({'x': 'ra', 'y': 'dec'})
        super().__init__(cat, parmap=parmap, *args, **kwargs)

    @require_validation(['x', 'y'])
    def _init_points(self):
        """
        In adition to regular points, this method inits SkyCoords
        """

        self._skypoints = SkyCoord(self._cat[self._parmap['x']], self._cat[self._parmap['y']], unit="deg")
        super()._init_points()


    def add_subregion(self, name, polytype = 'polygon', *args, **kwargs):
        """
        Shapely does not represent circular regions exactly, so this method
        allows for adding circular regions with a simple center-radius syntax
        Can create circles by passing polytype = 'circle', center = <float>, radius = <float>
        Center and radaius should be in units of degrees

        If no polytype is passed, assume it is a shapely.geometry.polygon
        """
        if polytype == 'polygon':
            super().add_subregion(name, *args, **kwargs)
        elif polytype == 'circle':
            try:
                center = kwargs['center']
                radius = kwargs['radius']
                self._add_circular_subregion(name, center, radius)
            except:
                logging.error("Expected a radius and a center, but couldn't find them.")

    def _add_circular_subregion(self, name, center, radius):
        """Internal function for initializing circular subregion"""
        centertype = type(center)
        if centertype == geometry.Point:
            input_center = SkyCoord(center.x, center.y, unit="deg")

        elif centertype == SkyCoord:
            input_center = center

        try:
            input_radius = radius.to(u.deg)
        except:
            logging.warning("No unit passed circle radius, assuming degrees")
            input_radius = radius*u.deg

        payload = {'polytype': 'circle', 'center': input_center, 'radius': input_radius}
        self._subregions.update({name: payload})

    @require_points
    def filter_by_subregion(self, name, *args, **kwargs):
        """
        Wrapper to deal with circular subregions implemented as center-radius pairs
        with SkyCoords.

        Otherwise just passes to the underlying Catalog2D method
        """
        try:
            region = self._subregions[name]
        except:
            logging.error("Subregion {} does not exist".format(name))
        try:
            mask = region['mask']
            return self._cat[mask]
        except:
            if region['polytype'] == 'circle':
                self._init_circular_subregion(name)
                return self._cat[self._subregions[name]['mask']]
            else:
                return super().filter_by_subregion(name, *args, **kwargs)

    def _init_circular_subregion(self, name):
        region = self._subregions[name]
        distances = region['center'].separation(self._skypoints)
        mask = distances <= region['radius']
        self._subregions[name]['mask'] = mask
