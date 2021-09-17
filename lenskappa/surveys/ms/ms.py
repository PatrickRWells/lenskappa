"""

The millenium simulation data is very different from a standard survey,
both in content and in the way its used.

possible todo - update from "survey" to dataset



"""

from hashlib import sha1
from lenskappa.region import CircularSkyRegion, SkyRegion
from shapely import geometry
from shapely.geometry import geo

from lenskappa.catalog import SkyCatalog2D
from lenskappa.simulations.simulation import Simulation
from lenskappa.params import QuantCatalogParam
from shapely import geometry
import pandas as pd


from shapely import geometry
import os
import re
import numpy as np
import logging
import astropy.units as u
from astropy.coordinates import SkyCoord


class millenium_simulation(Simulation):

    def __init__(self, *args, **kwargs):
        super().__init__("ms", *args, **kwargs)
        self._init_region()
        self._validate()

    def _init_region(self, *args, **kwargs):
        """
        The millenium simulation consists of 64 4x4 deg^2 fields
        each of which is subdividied into 16 1x1 deg^2 subfields.
        Because they are not actually on the sky, we can just
        initialize a single region and projected different catalogs
        onto it as needed
        """
        half_width = (2.0*u.degree).to(u.radian).value
        region = geometry.box(-half_width, -half_width, half_width, half_width) #in degrees
        center = geometry.Point(0, 0)
        self._region = SkyRegion(center, region, override=True, *args, **kwargs)
        self._init_subregions(*args, **kwargs)

    def _init_subregions(self, *args, **kwargs):
        width = (1.0*u.degree).to(u.radian).value
        min = (-1.5*u.degree).to(u.radian).value
        for x_i in range(4):
            for y_i in range(4):

                x_center = min + x_i*width
                y_center = min + y_i*width
                center = geometry.Point(x_center, y_center)
                x_min = x_center - width/2
                y_min = y_center - width/2
                x_max = x_center + width/2
                y_max = y_center + width/2
                subregion = geometry.box(x_min, y_min, x_max, y_max)
                key = "{}_{}".format(str(x_i), str(y_i))
                self._region.add_subregion(key, center, subregion, True, *args, **kwargs)


    def load_kappa_map(self, x: int, y: int, slice=36, filetype="binary"):
        """
        Loads the kappa map for a given field and slice.
        Expected size 4096/4096

        Params:
            x <int>: X-value of the field, should be between 0 and 7
            y <int>: Y-value of the field, should be between 0 and 7
            slice: Redshift slice number
        """

        map_directory = self._datamanager.get_file_location({'datatype': 'kappa_maps', 'slice': str(slice)})
        search_pattern = ".*8_{}_{}".format(str(x), str(y))
        files = [file for file in os.listdir(map_directory) if not file.startswith('.')]
        r = re.compile(search_pattern)
        matched_files = list(filter(r.match, files))
        if len(matched_files) > 1:
            logging.error("Too many kappa files found for index {} and {}".format(x, y))
            print(matched_files)
            return
        elif len(matched_files) == 0:
            logging.error("No kappa files found for index {} and {}".format(x,y))
            return

        kappa_data = self._load_kappa_file(os.path.join(map_directory, matched_files[0]))
        try:
            all_kappa = self._kappa_data
        except:
            self._kappa_data = {}

        self._kappa_data.update({"{}_{}".format(str(x),str(y)): kappa_data})

    def _load_kappa_file(self, file):
        """
        Kappa files are encodded as binary files with 4-byte floats
        Should always be 4096*4096
        """
        try:
            data = np.fromfile(file, dtype = np.float32)
        except:
            logging.error("Unable to load file {}".format(file))
            raise
        try:
            data = np.reshape(data, (4096,4096))
            return data
        except:
            logging.error("Unable to reshape kappa data into a 4096x4096 array")
            raise

    def load_catalogs_by_field(self, x, y, z_max = -1):
        """
        Loads catalogs for the field given by x,y
        """
        catalog_directory = self._datamanager.get_file_location({'datatype': 'catalogs', 'slice': 'global'})
        search_pattern = ".*8_{}_{}".format(str(x), str(y))
        files = [file for file in os.listdir(catalog_directory) if not file.startswith('.')]
        r = re.compile(search_pattern)
        matched_files = list(filter(r.match, files))
        if len(matched_files) > 16:
            logging.error("Too many catalog files found for index {} and {}".format(x, y))
            return
        elif len(matched_files) < 16:
            logging.error("Not enough catalog files found for index {} and {}".format(x,y))
            return
        cat = self._load_catalog_files(catalog_directory, matched_files, z_max=z_max)
        try:
            all_catalog = self._catalog_data
        except:
            self._catalog_data = {}
        key = "{}_{}".format(str(x), str(y))
        self._catalog_data.update({key: cat})

    
    def _load_catalog_files(self, directory, matched_files, z_max = -1):
        dfs = []
        for file in matched_files:
            subregion_key = re.search(r"\d_\d_\d_\d_\d", file).group()[-3:]
            try:
                df = pd.read_csv(os.path.join(directory, file), delimiter='\t')
                df['subregion'] = subregion_key
                dfs.append(df)
            except:
                logging.error("Unable to load file {}".format(file))
                raise
        
        combined = pd.concat(dfs, ignore_index=True)
        ra_par = QuantCatalogParam("pos_0[rad]", 'ra', u.radian)
        dec_par = QuantCatalogParam("pos_1[rad]", "dec", u.radian)
        pars = [ra_par, dec_par]
        if z_max > 0:
            combined.drop(combined[combined['z_spec'] > z_max].index, inplace=True)
        return SkyCatalog2D(combined, params=pars)


    def compute_weights(self, x, y, z_max, aperture = 120*u.arcsec, *args, **kwargs):
        """
        Compute weights for a given field
        """
        x_grid, y_grid = self._generate_grid(aperture, *args, **kwargs)
        self.load_catalogs_by_field(x, y, z_max = z_max)
        key = "{}_{}".format(str(x),str(y))
        catalog = self._catalog_data[key]
        for i_x in x_grid:
            for i_y in y_grid:
                center = millenium_simulation.get_position_from_index(i_x, i_y)
                newcat = self._filter_catalog(catalog, center, aperture)
                print("{}, {}: {}".format(str(i_x), str(i_y), len(newcat)))
            


    def _generate_grid(self, aperture = 120*u.arcsec, *args, **kwargs):
        #First, find the corners of the tiling region.
        #Since we don't allow tiles to overlap with the edge of the field.
        min_pos = -2.0*u.degree + aperture
        max_pos = 2.0*u.degree - aperture
        bl_corner = millenium_simulation.get_index_from_position(min_pos, min_pos)
        tr_corner = millenium_simulation.get_index_from_position(max_pos, max_pos)
        #Since there's rounding involved above, check to make sure the tiles don't
        #Overlap with the edge of the field.
        min_pos_x, min_pos_y = millenium_simulation.get_position_from_index(*bl_corner)
        max_pos_x, max_pos_y = millenium_simulation.get_position_from_index(*tr_corner)

        min_vals = -2.0*u.degree
        max_vals = 2.0*u.degree
        
        x_diff = min_pos_x - min_vals
        y_diff = min_pos_y - min_vals
        x_index = bl_corner[0]
        y_index = bl_corner[1]
        if x_diff < aperture:
            x_index += 1
        if y_diff < aperture:
            y_index += 1
        bl_corner = (x_index, y_index)

        x_diff = max_vals - max_pos_x
        y_diff = max_vals - max_pos_y
        x_index = tr_corner[0]
        y_index = tr_corner[1]
        if x_diff < aperture:
            x_index -= 1
        if y_diff < aperture:
            y_index -= 1
        tr_corner = (x_index, y_index)

        min_pos_x, min_pos_y = millenium_simulation.get_position_from_index(*bl_corner)
        max_pos_x, max_pos_y = millenium_simulation.get_position_from_index(*tr_corner)

        x_pos = min_pos_x
        x_grid = []
        while x_pos < max_pos_x:
            i_x, i_y = millenium_simulation.get_index_from_position(x_pos, min_pos_y)
            x_grid.append(i_x)
            x_pos += 2*aperture
        y_pos = min_pos_y
        y_grid = []
        while y_pos < max_pos_y:
            i_x, i_y = millenium_simulation.get_index_from_position(min_pos_x, y_pos)
            y_grid.append(i_y)
            y_pos += 2*aperture
        
        return x_grid, y_grid

    def _filter_catalog(self, cat, center, aperture):
        cat.get_points()
        poly_r = aperture.to(u.radian).value
        poly = geometry.Point(center[0].value, center[1].value).buffer(poly_r)
        subregion_overlaps = self._region.get_subregion_intersections(poly)
        region = CircularSkyRegion(SkyCoord(*center), aperture)
        newcat = cat.filter_by_columns([{'subregion': subregion_overlaps}])
        final_cat = newcat.get_objects_in_region(region)
        return final_cat
        
    @classmethod
    def get_position_from_index(cls, x, y):
        """
        Returns an angular position (in radians) based on a given x, y index
        Where x,y are in the range [0, 4096]. This matches with the
        kappa maps
        The position returned is with reference to the center of the field,
        so negative values are possible
        """
        l_field = (4.0*u.degree).to(u.radian) #each field is 4 deg x 4 deg
        n_pix = 4096.0
        l_pix = l_field/n_pix
        pos_x =  -0.5*l_field + (x+0.5) * l_pix
        pos_y =  -0.5*l_field + (y+0.5) * l_pix
        return pos_x,pos_y
    
    @classmethod
    def get_index_from_position(cls, pos_x, pos_y):
        """
        Returns the index of the nearest kappa point given an angular position.
        
        """
        try:
            pos_x_rad = pos_x.to(u.radian)
            pos_y_rad = pos_y.to(u.radian)
        except:
            logging.error("Need angular distances to get kappa indices!")
            raise

        l_field = (4.0*u.degree).to(u.radian) #each field is 4 deg x 4 deg
        n_pix = 4096.0
        l_pix = l_field/n_pix

        x_pix = ( (pos_x_rad + 0.5*l_field) / l_pix) - 0.5
        y_pix = ( (pos_y_rad + 0.5*l_field) / l_pix) - 0.5
        return int(round(x_pix.value)), int(round(y_pix.value))







if __name__ == "__main__":
    ms = millenium_simulation()
    ms.compute_weights(1, 1, z_max = 1.523)