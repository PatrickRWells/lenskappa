from astropy.coordinates.sky_coordinate import SkyCoord, SkyOffsetFrame
from numpy.lib.polynomial import poly
import pandas as pd
import numpy as np
from astropy import wcs, coordinates, units
import astropy.units as u
from shapely import geometry
from copy import copy, deepcopy
import logging
import shapely
from shapely.geometry import geo
import toml
from abc import ABCMeta, abstractmethod
from lenskappa.region import CircularSkyRegion, SkyRegion

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
                logging.warning("Catalog parameters have not been validated!")
                return None
            isvalid = [checkvalid[name] for name in params]
            if np.all(isvalid):
                return fn(self, *args, **kwargs)
            else:
                logging.error("Catalog parameter(s) {} are not valid".format([par for index, par in enumerate(pars) if not isvalid[index]]))
        return wrapper
    return outer

class Catalog(metaclass=ABCMeta):

    def __init__(self, cat, base_parmap = {}, partypes = {}, *args, **kwargs):
        self._cat = cat
        self._needed_parameter_types = partypes
        self._parmap = base_parmap
        try:
            new_parmap = kwargs['parmap']
            self.load_params(new_parmap)
        except:
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

    def get_dataframe(self):
        return self._cat

    def get_parmap(self):
        return self._parmap

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
    @classmethod
    def from_dataframe(cls, cat, *args, **kwargs):
        return cls(cat, *args, **kwargs)
    
    def apply_catalog_mask(self, mask):
        df = self._cat.loc[self._cat[mask].index]
        return self.from_dataframe(df, parmap=self._parmap)

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

        for par, parmap in input_map.items():
            self._parmap.update({par: parmap})
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
    def __init__(self, cat, base_parmap = {}, *args, **kwargs):
        needed_parameter_types = {'x': float, 'y': float}
        super().__init__(cat, base_parmap, partypes=needed_parameter_types, *args, **kwargs)
        self._points_initialized = False
        self._subregions = {}
        self._unions = {}


    def get_points(self):
        if not self._points_initialized:
            self._init_points()
        return self._points

    def add_subregion(self, name, reg, *args, **kwargs):
        """
        Adds a subregion to the catalog
        This is mostly just a convinient way of subdividing catalog data

        params:
            name <str>: Name for the subregion
            polygon <shapely.geometry.Polygon>: Polygon defining the region

        """
        regtype = type(reg)
        if regtype in [SkyRegion, CircularSkyRegion]:
            self._subregions.update({name: {'region': reg} } )
        else:
            logging.error("Tried to add a subregion, but the region was not a polygon")    


    @require_points
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
        
        try:
            mask = self._subregions[name]['mask']
        except:
            mask = self._init_subregion[name]

        return self.from_dataframe(cat=self._cat[mask], parmap=self._parmap)
    
    def add_column_bound(self, column, min = False, max = False):
        """
        Remove objects outside of a particular bound, based on value
        """
        try:
            colname = self._parmap[column]
        except:
            colname = column
        try:
            col = self._cat[colname]
        except:
            logging.error("Tried to filter catalog by column {}, but the column doesn't exist")

        if not (min or max) and min !=0 and max !=0:
            logging.error("You must provide a minimum or maximum to place bounds on the catalog!")
            return
        
        mask = np.array([True]*len(self._cat))
        if min or min == 0:
            mask = mask & (col >= min)
        if max or max == 0:
            mask = mask & (col <= max)
        
        self._cat.drop(self._cat[~mask].index, inplace=True)




    def filter_by_columns(self, conditions, *args, **kwargs):
        """
        Filters the catalog by certain columns.
        Shold be key-value pairs, where the key is the column
        And the value is a list of acceptable values
        """
        final_mask = np.array([False]*len(self))
        for data in conditions:
            sub_mask = np.array([True]*len(self))
            for column, keys in data.items():
                sub_mask = sub_mask & self._cat[column].isin(keys)

            final_mask = final_mask | sub_mask
        
        return self.apply_catalog_mask(final_mask)


    @require_points
    def filter_by_subregions(self, subregions, *args, **kwargs):

        masks = []
        for region_id in subregions:
            try:
                subregion = self._subregions[region_id]
            except:
                logging.error("Subregion {} does not exist. Skipping...")
                continue
            try:
                masks.append(subregion['mask'])
            except:
                masks.append(self._init_subregion(region_id))
        if masks:
            final_mask = np.any(masks, axis=0)
            return self.apply_catalog_mask(mask = final_mask)
        else:
            logging.error("None of the subregions were found.")


    @require_validation(params=['x', 'y'])
    def _init_points(self):
        """
        Initialize a set of shapely points and store them in the catalog
        TODO: Stress-test with millions of points
        """

        logging.info("Initializing catalog points")
        x = self._cat[self._parmap['x']]
        y = self._cat[self._parmap['y']]
        #This line takes quite a while for large catalogs
        point_objs = [geometry.Point(xy) for xy in list(zip(x,y))]
        self._points = pd.Series(point_objs)
        self._points_initialized = True

    @require_points
    def _init_subregion(self, name):
        """
        Initializes subregion by finding which points fall inside of it
        Results are stored as a mask that can be applied to _cat
        Method assumes that _init_points has already been run
        """
        poly = self._subregions[name]
        #Type validation of polygons is handled when the subregion is created
        mask = [point.within(poly) for point in self._points]
        initial_len = len(self._cat)
        final_len = len(self._cat[mask])
        logging.info("Region {} masks out {} objects".format(name, initial_len-final_len))
        self._subregions[name]['mask'] = mask
        return mask

    @require_points
    def get_objects_in_region(self, region):
        #Note, this is a very simple implementationt that will be inefficient for very large regions
        #It is recommended that you override this for larger catalogs 
        mask = np.array([region.contains(point) for point in self._points])
        return self[mask]
        


class SkyCatalog2D(Catalog2D):

    def __init__(self,cat, *args, **kwargs):
        base_parmap = ({'x': 'ra', 'y': 'dec'})
        super().__init__(cat, base_parmap=base_parmap, *args, **kwargs)

    @require_validation(['x', 'y'])
    def _init_points(self):
        """
        In adition to regular points, this method inits SkyCoords
        """

        #this is very slow. Need to work out a different solution
        super()._init_points()

    def get_coords(self):
        try:
            return self._skypoints
        except:
            self._skypoints = SkyCoord(self._cat[self._parmap['x']], self._cat[self._parmap['y']], unit="deg")
            return self._skypoints

    
    def add_subregion(self, name, region, *args, **kwargs):
        """
        Shapely does not represent circular regions exactly, so this method
        allows for adding circular regions with a simple center-radius syntax
        Can create circles by passing polytype = 'circle', center = <float>, radius = <float>
        Center and radaius should be in units of degrees

        If no polytype is passed, assume it is a shapely.geometry.polygon
        """
        self._polytype = region.get_type()
        super().add_subregion(name, region)

    def rotate(self, original, new, *args, **kwargs):
        """
        Rotates a catalog on the sky
        """
        if len(self) > 30000:
            logging.error("Tried to rotate this catalog, but its too big!")
            return
        coords = self.get_coords()
        separations = original.separation(coords)
        pas = original.position_angle(coords)
        new_coords = new.directional_offset_by(pas, separations)

        new_df = deepcopy(self._cat)
        new_df[self._parmap['x']] = new_coords.ra.degree
        new_df[self._parmap['y']] = new_coords.dec.degree
        new_catalog = self.from_dataframe(new_df, parmap=self._parmap)
        new_catalog._skypoints = new_coords
        return new_catalog

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

    def _init_circular_subregion(self, name):
        region = self._subregions[name]
        distances = region['center'].separation(self._skypoints)
        mask = distances <= region['radius']
        self._subregions[name]['mask'] = mask

    def get_objects_in_region(self, region, *args, **kwargs):
        if type(region) == CircularSkyRegion:
            center, radius = region.skycoord

            coords = self.get_coords()
            distances = coords.separation(center)
            self._cat['dist'] = distances
            mask = (distances <= radius)
            self._parmap.update({'r': 'dist'})
            item = self.apply_catalog_mask(mask)
            return item
        else:
            return super().get_objects_in_region(region)
