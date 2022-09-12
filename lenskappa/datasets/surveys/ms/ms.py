"""

The millenium simulation data is very different from a standard survey,
both in content and in the way its used.

possible todo - update from "survey" to dataset



"""

from genericpath import exists
from importlib.resources import path
from math import gamma
from nis import match
from turtle import st
from lenskappa.spatial import CircularSkyRegion, SkyRegion
from shapely import geometry

from lenskappa.catalog import SkyCatalog2D
from lenskappa.datasets.simulations.simulation import Simulation
from lenskappa.catalog import QuantCatalogParam, SingleValueParam
from shapely import geometry
import pandas as pd


from shapely import geometry
import os
import re
import numpy as np
import logging
import astropy.units as u
from astropy.coordinates import SkyCoord
import pathlib
import multiprocessing
from functools import partial

class millenium_simulation(Simulation):

    def __init__(self, *args, **kwargs):

        """
        DataSet class used for the millenium simulation.
        
        """
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
        width = (4.0*u.degree).to(u.radian).value
        region = geometry.box(0, 0, width, width) #in degrees
        center = geometry.Point(2, 2)
        self._region = SkyRegion(center, region, override=True, *args, **kwargs)
        self._init_subregions(*args, **kwargs)

    def _init_subregions(self, *args, **kwargs):
        """
        Each millenium simulation field is broken up into 16 1deg^2 fields. 
        While there's nothing special about them, they are useful for simplifying
        object finding.        
        """

        width = (1.0*u.degree).to(u.radian).value
        min = (0.5*u.degree).to(u.radian).value
        #Locations in the MS are given in radians for reasons unknown.
        #Instead of using the standared -2 deg -> 2-deg, I reindex to
        #0deg -> 4deg

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
        Loads the binary kappa maps for a given field and slice.
        Expected size 4096/4096. 
        Expects the kappa maps to be added in the datamanager

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
        Loads a given kappa file into a 2D numpy array.
        Kappa files are encodded as binary files with 4-byte floats
        Should always be 4096 x 4096
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
    
    def load_gamma_maps(self, x: int, y: int, slice=36, filetype="binary"):
        map_directory = self._datamanager.get_file_location({'datatype': 'gamma_maps', 'slice': str(slice)})
        search_pattern = ".*8_{}_{}".format(str(x), str(y))
        files = [file for file in os.listdir(map_directory) if not file.startswith('.')]
        r = re.compile(search_pattern)
        matched_files = list(filter(r.match, files))
        if len(matched_files) > 1:
            logging.error("Too many kappa files found for index {} and {}".format(x, y))
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




    def load_catalogs_by_field(self, x, y, z_s = -1, params = []):
        """
        Loads catalogs for the field given by x,y
        Catalogs are expected to be added to the datamanger
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
        cat = self._load_catalog_files(catalog_directory, matched_files, z_max=z_s, params = params)
        z_s_par = SingleValueParam("z_s", z_s)
        cat.add_param(z_s_par)
        self._catalog = cat

    def attach_dist_array(self, dist_array, target_param):
        self._catalog.attach_dist_array(dist_array, target_param)


    def _load_catalog_files(self, directory, matched_files, z_max = -1, params = []):
        """
        The data for each of the 64 fields is broken into 1x1 deg fields
        And each of those has their own catalog. Here, we load them
        all into a single dataframe.
        """

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
        z_par = QuantCatalogParam("z_spec", "z_gal")

        pars = [ra_par, dec_par, z_par] + params
        #Re-index positions from -2deg -> 2deg to 0deg->4deg
        combined['pos_0[rad]'] += (2.0*u.degree).to(u.radian).value
        combined['pos_1[rad]'] += (2.0*u.degree).to(u.radian).value

        if z_max > 0: #Drop objects above the maximum redshift
            combined.drop(combined[combined['z_spec'] > z_max].index, inplace=True)
        return SkyCatalog2D(combined, params=pars)


    def _generate_grid(self, aperture = 120*u.arcsec, overlap = 1, *args, **kwargs):
        """
        Generates a grid of locations to compute weighted number counts on.
        Here, the locations are determined by the grid points defined in the
        millenium simulation.

        Params:

        aperture: Size of the aperture to consider. Should be an astropy quantity
        overlap: If 1, adjacent tiles do not overlap (centers have spacing of
            2*aperture). If above 1, tile spacing = 2*(aperture/overlap).
        """


        #First, find the corners of the tiling region.
        #Since we don't allow tiles to overlap with the edge of the field.
        min_pos = 0.0*u.degree + aperture
        max_pos = 4.0*u.degree - aperture
        bl_corner = millenium_simulation.get_index_from_position(min_pos, min_pos)
        tr_corner = millenium_simulation.get_index_from_position(max_pos, max_pos)
        #Since there's rounding involved above, check to make sure the tiles don't
        #Overlap with the edge of the field.
        min_pos_x, min_pos_y = millenium_simulation.get_position_from_index(*bl_corner)
        max_pos_x, max_pos_y = millenium_simulation.get_position_from_index(*tr_corner)

        min_vals = 0.0*u.degree
        max_vals = 4.0*u.degree
        pix_distance = 4.0*u.deg/4096.0

        x_diff = min_pos_x - min_vals
        y_diff = min_pos_y - min_vals
        x_index = bl_corner[0]
        y_index = bl_corner[1]

        #Make sure we're fully within the field
        if x_diff < aperture:
            x_index += 1
        if y_diff < aperture:
            y_index += 1
        bl_corner = (x_index, y_index)

        x_diff = max_vals - max_pos_x
        y_diff = max_vals - max_pos_y
        x_index = tr_corner[0]
        y_index = tr_corner[1]

        #Make sure we're fully within the field.
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
            if overlap == 'all':
                x_pos += pix_distance
            else:
                x_pos += 2* (aperture/overlap)
        y_pos = min_pos_y
        y_grid = []
        while y_pos < max_pos_y:
            i_x, i_y = millenium_simulation.get_index_from_position(min_pos_x, y_pos)
            y_grid.append(i_y)
            if overlap == 'all':
                y_pos += pix_distance
            else:
                y_pos += 2* (aperture/overlap)
        
        return x_grid, y_grid

    def get_ciruclar_tile(self, aperture, *args, **kwargs):
        """
        Generator that returns ciruclar tiles, one at a time.
        
        Params:

        Aperture: Size of the tiles
        
        """

        x_grid, y_grid = self._generate_grid(aperture, *args, **kwargs)
        #Assumption is that the aperture isn't changing during a run
        #maybe should change that

        for x_i in x_grid:
            for y_i in y_grid:
                center = millenium_simulation.get_position_from_index(x_i, y_i)
                region = CircularSkyRegion(SkyCoord(*center), aperture)
                yield region

    def has_samples(self, *args, **kwargs):
        return self._catalog.has_samples(*args, **kwargs)


    def get_objects_in_region(self, region):
        """
        Gets the objects in a particular region.
        """
        self._catalog.get_points()
        poly = region.get_polygon(unit=u.radian)
        subregion_overlaps = self._region.get_subregion_intersections(poly)
        newcat = self._catalog.filter_by_columns([{'subregion': subregion_overlaps}])
        final_cat = newcat.get_objects_in_region(region)
        return final_cat

        
    @classmethod
    def get_position_from_index(cls, x, y):
        """
        Returns an angular position (in radians) based on a given x, y index
        Where x,y are in the range [0, 4096]. This matches with the
        gridpoints defined by the millenium simulation.

        The position returned is with reference to the center of the field,
        so negative values are possible
        """
        l_field = (4.0*u.degree).to(u.radian) #each field is 4 deg x 4 deg
        n_pix = 4096.0
        l_pix = l_field/n_pix
        pos_x =  (x+0.5) * l_pix
        pos_y =  (y+0.5) * l_pix
        return pos_x,pos_y
    
    @classmethod
    def get_index_from_position(cls, pos_x, pos_y):
        """
        Returns the index of the nearest grid point given an angular position.
        
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

        x_pix = pos_x/l_pix - 0.5
        y_pix = pos_y/l_pix -0.5
        return int(round(x_pix.value)), int(round(y_pix.value))
    
def extract_maps(input_directory: pathlib.Path, output_directory: pathlib.Path, threads = 1, override = False) -> None:
    """
    Converts original .Phi files from Stefan Hilbert into kappa and gamma maps.
    Input: directory containing original maps
    Output: Outuput directory. Will create two subdirectories, one for kappa and one for gamma
    """
    indices = list(range(8))
    files = [p.name for p in input_directory.glob('*.Phi')]
    #print(files)
    plane_re = lambda x: re.search(r'plane_\d+', x)
    matches = [plane_re(f) for f in files]
    planes = [int(str(m.group(0)).split("_")[1]) for m in matches]
    if not all(p == planes[0] for p in planes):
        print("Warning: Not all the .Phi files come from the same plane!")
        if not override:
            print("Set override = True to continue anyway")
            return
        planes = np.unique(planes)
        print("Subdirectories will be created for each plane")
        for p_ in planes:
            plane_output = output_directory / f"plane_{p_}"
            plane_output.mkdir(exist_ok=False)
            matching_files = filter(files, lambda x: f"plane_{p_}" in x)
            for f in matching_files:
                os.symlink(input_directory / f, plane_output / f)
            millenium_simulation.extract_maps(plane_output, plane_output)
            for f in matching_files:
                pathlib.unlink(plane_output / f)
            
            print("Done!")
            return
    if(len(files)) != 64:
        print("Warning: Missing files for this plane!")

    kappa_output = output_directory / "kappa"
    gamma_output = output_directory / "gamma"
    kappa_output.mkdir(exist_ok = True)
    gamma_output.mkdir(exist_ok = True)
    print("Extracting convergence and shear information, this will take some time")
    if threads == 1:
        print("If you like, you can speed this up by chanigng the \"threads\" parameter")
    if len(files) < threads:
        threads = len(files)
    
    with multiprocessing.Pool(threads) as pool:
        f = partial(_extract_single_file, input_dir = input_directory, kappa_output = kappa_output, gamma_output = gamma_output)
        pool.map(f, files)


def _extract_single_file(f, input_dir, kappa_output, gamma_output):
    print(f"Extracting file {f}")
    data = np.fromfile(input_dir / f, np.float32)
    kappas = np.ndarray((4096, 4096), np.float32)        
    gamma1 = np.ndarray((4096, 4096), np.float32)        
    gamma2 = np.ndarray((4096, 4096), np.float32)
    for i in range(4096):
        if i%1000 == 0:
            print(f"File {f}: extracted {i}/4096 points")
        for j in range(4096):
            start = 4*4096*i + 4*j
            m = data[start: start+4]
            k_, g1_, g2_ = _extract_single_matrix(m)
            kappas[i,j]  = k_
            gamma1[i,j] = g1_
            gamma2[i,j] = g2_
    fname = pathlib.Path(f).stem
    kappa_out = kappa_output / ".".join([fname, 'kappa'])

    g1_out = gamma_output / ".".join([fname, 'gamma_1'])
    g2_out = gamma_output / ".".join([fname, 'gamma_2'])
    kappas.tofile(kappa_out)
    gamma1.tofile(g1_out)
    gamma2.tofile(g2_out)
    print(f"Finished extracting file{f}")

def _extract_single_matrix(mat: list) -> tuple:
    a,b,c,d = mat
    tr = a + c
    btr = b - d
    s = np.sign(tr)
    hyp = np.sqrt(tr*tr + btr*btr)
    kappa = kappa = 1.0 - 0.5*s*hyp

    dn = -0.5*s / hyp
    g1 = dn * (a*a - b*b + c*c - d*d)
    g2 = dn * 2 * (a*b+c*d)
    return kappa, g1, g2
