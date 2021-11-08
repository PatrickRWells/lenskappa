from lenskappa.params import QuantCatalogParam
import pandas as pd
import numpy as np
import astropy.units as u
import logging

from astropy.coordinates.sky_coordinate import SkyCoord
from shapely import geometry
from abc import ABC, abstractmethod

from lenskappa.region import CircularSkyRegion, SkyRegion
from lenskappa.utils.decorators import require_points, require_validation
from lenskappa.params import CatalogParam


class Catalog(ABC):

    def __init__(self, cat, base_parmap = {}, partypes = {}, *args, **kwargs):
        self._cat = cat
        self._needed_parameter_types = partypes
        self._base_parmap = base_parmap
        try:
    
            new_parmap = kwargs['params']
            self.load_params(new_parmap)
        except:
            try:
                self._parmap = kwargs['parmap']
            except:
                pass
            self._validate_parmap()

    def __getitem__(self, key):
        """
        Allows for masking as with a usual dataframe
        """
        if key in self._parmap.keys():
            return self._parmap[key].get_values(self._cat)
        elif key in self._inverse_map.keys():
            return self._parmap[self._inverse_map[key]].get_values(self._cat)
        elif key in self._base_parmap.keys():
            return self[self._base_parmap[key]]
        else:
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

    def get_colname(self, column):
        if column in self._cat.columns:
            return column
        try:
            return self._parmap[column].col
        except:
            logging.warning("Unable to find column {}".format(column))
            raise



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
            raise

    @classmethod
    def from_dataframe(cls, cat, *args, **kwargs):
        return cls(cat, *args, **kwargs)

    def apply_boolean_mask(self, mask):
        df = self._cat.loc[self._cat[mask].index]
        try:
            input_points = self._points[mask]
            return self.from_dataframe(df, parmap=self._parmap, points = input_points)
        except:
            return self.from_dataframe(df, parmap=self._parmap)


    def replace_values_by_mask(self, mask, column, value):
        """
        Returns a new catalog where values that fall outside the mask
        Are replaced with the given value.
        """
        df = self._cat.copy()
        colname = self.get_colname(column)
        df.loc[~mask, colname] = value
        try:
            input_points = self._points[mask]
            return self.from_dataframe(df, parmap=self._parmap, points = input_points)
        except:
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
            parmap = self._base_parmap
        except:
            self._base_parmap = {}
        inputed_pars = [par.standard for par in input_map]
        for par, parmap in self._base_parmap.items():
            if parmap in inputed_pars:
                pass
            else:
                logging.warning("No rename found for parameter {}, defaulting to {}".format(par, parmap))
        self._parmap = {}
        for par in input_map:
            self._parmap.update({par.standard: par})

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
            try:
                single_par = self._parmap[par]
            except:
                base_par = self._base_parmap[par]
                if base_par not in self._parmap.keys():
                    try:
                        col = self._cat[base_par]
                        single_par = CatalogParam(base_par, base_par)
                    except:
                        logging.error("Unable to find column {} needed for this catalog type!".format(single_par))
                        raise
                else:
                    single_par = self._parmap[base_par]

            try:
                col = single_par.get_values(self._cat)
            except:
                logging.warning("Unable to find column {} in data".format(single_par.col))
                self._valid.update({par: False})
                continue

            try:
                coldata = col.astype(partype)
                self._valid.update({par : True})
            except:
                logging.error("Unable to cast column {} as type {}".format(single_par.col, partype))
                self._valid.update({par : False})

            self._inverse_map = {p.col: p.standard for p in self._parmap.values()}



    def add_param(self, par):
        name = par.standard
        if name not in self._parmap.keys():
            self._parmap.update({name: par})

    def add_params(self, pars):
        for par in pars:
            self.add_param(par)

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
        self._subregions = {}
        self._unions = {}


    def get_points(self, *args, **kwargs):
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
        try:
            input_points = self._points[mask]
            return self.from_dataframe(cat=self._cat[mask], parmap=self._parmap, points = input_points)
        except:
            return self.from_dataframe(cat=self._cat[mask], parmap=self._parmap)

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
            return self.apply_boolean_mask(mask = final_mask)
        else:
            logging.error("None of the subregions were found.")


    @require_validation(params=['x', 'y'])
    def _init_points(self):
        """
        Initialize a set of shapely points and store them in the catalog
        TODO: Stress-test with millions of points
        """

        x = self['x'].value
        y = self['y'].value
        #This line takes quite a while for large catalogs
        point_objs = [geometry.Point(*xy) for xy in list(zip(x,y))]
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

    def get_skypoints(self):

        try:
            return self._skypoints
        except:
            self._skypoints = SkyCoord(self['ra'], self['dec'])
            self._skypoints_initialized = True
            return self._skypoints
    
    def _get_points(self):
        try:
            p = self._points
            return p
        except AttributeError:
            self._init_points()
            self._points_initialized
            return self._points

    def get_points(self, point_type='shapely'):
        if point_type == 'shapely':
            return self._get_points()
        else:
            return self._skypoints()

    def get_distances(self, center: SkyCoord, unit=u.arcsec, *args, **kwargs):
        coords = self.get_skypoints()
        distances = center.separation(coords).to(unit).value
        self._cat['distance'] = distances
        dist_par = QuantCatalogParam('distance', 'r', unit)
        self._parmap.update({dist_par.standard: dist_par})


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
        coords = self.get_skypoints()
        separations = original.separation(coords)
        pas = original.position_angle(coords)
        new_coords = new.directional_offset_by(pas, separations)

        new_df = self._cat.copy()
        new_df[self._parmap['ra']] = new_coords.ra.degree
        new_df[self._parmap['dec']] = new_coords.dec.degree
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
            self.get_distances(center)
            mask = self['r'] <= radius
            item = self.apply_boolean_mask(mask)
            return item
        else:
            return super().get_objects_in_region(region)

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

        return self.apply_boolean_mask(final_mask)
