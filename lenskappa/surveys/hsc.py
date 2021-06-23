from numpy.lib.polynomial import poly
import shapely
from lenskappa import region
from lenskappa.surveys.survey import Survey
from lenskappa.region import SkyRegion
import re
from astropy.coordinates import SkyCoord
from astropy import wcs
import astropy.units as u
from shapely import geometry
import pandas as pd
import numpy as np
import math
import os
import time
import functools
from copy import copy
import logging
import toml
import operator

class HSCSurvey(Survey):
    def __init__(self):
        pass

    @classmethod
    def read_hsc_tracts(cls, tractfile):
        tracts = HSCSurvey._parse_tractfile(tractfile)
        regions = HSCSurvey._parse_tractdata(tracts)
        return regions
    
    @staticmethod
    def _parse_tractfile(tractfile):
        tracts = {}
        with open(tractfile) as tf:
            for line in tf:
                if line.startswith('*'):
                    continue

                nums = re.findall(r"[-+]?\d*\.\d+|\d+", line)
                tract_num = int(nums[0])
                patch = re.search("Patch", line)
                if patch is None:
                    if tract_num not in tracts.keys():
                        tracts.update({tract_num: {'corners': [], 'type': 'tract', 'subregions': {} } } )
                    if 'Corner' in line:
                        tracts[tract_num]['corners'].append((float(nums[-2]), float(nums[-1])))
                    elif 'Center' in line:
                        tracts[tract_num].update({'center': (float(nums[-2]), float(nums[-1]))})

                else:
                    patch = re.findall(r'\d,\d', line)
                    patch_val = tuple(map(int, patch[0].split(',')))
                    if patch_val not in tracts[tract_num]['subregions'].keys():
                        tracts[tract_num]['subregions'].update({patch_val: {'corners': [], 'type': 'patch'}})
                    if 'Corner' in line:
                        tracts[tract_num]['subregions'][patch_val]['corners'].append((float(nums[-2]), float(nums[-1])))
                    elif 'Center' in line:
                        tracts[tract_num]['subregions'][patch_val].update({'center': (float(nums[-2]), float(nums[-1]))})
        return tracts
    
    @staticmethod
    def _parse_tractdata(tractdata):
        output = {}
        for name, tract in tractdata.items():
            corners = tract['corners']
            center = tract['center']
            polygon = HSCSurvey._parse_polygon_corners(center, corners)
            region_obj = SkyRegion(center, polygon)
            for patchname, patch in tract['subregions'].items():
                patch_corners = patch['corners']
                patch_center = patch['center']
                patch_polygon = HSCSurvey._parse_polygon_corners(patch_center, patch_corners)
                added = region_obj.add_subregion(HSCSurvey._patch_tuple_to_int(patchname), center, patch_polygon, override = True)
                if not added:
                    logging.warning("Failed to add subregion to tract {}".format(name))
            output.update({name: region_obj})
        return output
    
    @staticmethod
    def _patch_tuple_to_int(patch_tuple):
        if patch_tuple[0] == 0:
            return patch_tuple[1]
        else:
            return 100*patch_tuple[0] + patch_tuple[1]



    @staticmethod
    def _parse_polygon_corners(center, points):
        sorted_coords = sorted(points, key=lambda coord: (-135 - math.degrees(math.atan2(*tuple(map(operator.sub, coord, center))[::-1]))) % 360)

        #This ensures the points are ordered counter-clockwsie, to avoid twisted polygons
        #Shoutout to StackOverflow
        polygon = geometry.Polygon(sorted_coords)
        return polygon



class hsc_region_parser:
    def __init__(self, f):
        self._f = f
    def __call__(self, instance, subregion_name, *args, **kwargs):
        print("Changing subregion name")
        new_subregion_name = self._parse_subregion(subregion_name)
        self._f(instance, new_subregion_name, *args, **kwargs)
    def _parse_subregion(self, subregion_name):
        if type(subregion_name) == tuple:
            return self._parse_patch_tuple_int(subregion_name)
        elif type(subregion_name) == int:
            return subregion_name
        else:
            logging.warning("Error: Tried to parse HSC region name. Expected a tupe or an int, got {}".format(subregion_name))
            return subregion_name

    def __get__(self, instance, owner):
        return functools.partial(self, instance) 


def hsc_config(f):
    def wrapper(self, *args, **kwargs):
        config = toml.load("config/hsc.toml")
        kwargs.update({"config": config})
        return f(self, *args, **kwargs)
    return wrapper
        



class hsc_catalog:
    def __init__(self, file):
        self._data = pd.read_csv(file)
        print("Splitting data...")
        self.subsets = {id: self._data[self._data.tract == id] for id in self._data.tract.unique()}
        self._patch_data = {}

    
    def get_objects_within_patches(self, tracts):
        frames = []
        for id, data in tracts.items():
            patch_ints = [self._patch_tuple_to_int(key) for key in data['subregions'].keys()]
            for key in patch_ints:
                if key not in self._patch_data.keys():
                    mask = (self.subsets[id].patch).isin(patch_ints)
                    self._patch_data.update({key: self.subsets[id][mask]})
                frames.append(self._patch_data[key])
        return pd.concat(frames)

    
    def apply_mask_to_objects(self, patches, mask):
        objects = self.get_objects_within_patches(patches)
        return objects

class field:
    def __init__(self, top_region, sub_regions, catalog = None):
        self._catalog = catalog
        self._region = top_region
        self._subregions = {key: region(value) for key,value in sub_regions.items()}
        self._tiled = False
    def get_inside(self, center, aperture):
        #Given a circular aperture, determine which subregions it intersects with
        subregions = {}
        aperture_deg = aperture.to(u.deg).value
        circle_center = geometry.point.Point(center.ra.degree, center.dec.degree)
        aperture_poly = circle_center.buffer(aperture_deg)
        for id, region in self._subregions.items():
            inside = region.get_inside(aperture_poly)
            if inside or inside is None:
                subregions.update({id: inside})
        return subregions
    
    def get_weight_ratios(self, x, y, ext_catalog = None, ext_catalog_center = None, weights = None, hsc_params = {}, ext_params = {}):
        regmask = self.get_region_mask(x,y)
        center = SkyCoord(self._x_tiles[x], self._y_tiles[y], unit="deg")
        subregions = self.get_inside(center, self._aperture_size)
        local_cat = self._catalog.get_objects_within_patches(subregions)
        
        cat = self.apply_region_mask(center, self._aperture_size, regmask, local_cat, get_dist=True)

        new_ext_cat = ext_catalog.copy(deep=True)
        dra, ddec = center.ra.degree - ext_catalog_center.ra.degree, center.dec.degree - ext_catalog_center.dec.degree
        
        new_ext_cat.ra += dra
        new_ext_cat.dec += ddec

        ext_cat_masked = self.apply_region_mask(center, self._aperture_size, regmask, new_ext_cat)

        field_weights = {key: weight.compute_weight(ext_cat_masked, ext_params) for key, weight in weights.items()}
        control_weights = {key: weight.compute_weight(cat, hsc_params) for key, weight in weights.items()}

    def get_region_mask(self, x, y):
        if not self._tiled:
            return
        center = SkyCoord(self._x_tiles[x], self._y_tiles[y], unit="deg")        
        bright_mask_objects =  self.get_bright_mask(center, self._aperture_size)
        return bright_mask_objects


    def apply_region_mask(self, center, aperture, regmask, catalog, get_dist = False):
        object_coords = SkyCoord(catalog.ra, catalog.dec, unit="deg")
        seps = object_coords.separation(center)
        if get_dist:
            catalog['dist'] = seps.to(u.arcsec)

        in_aperture_mask = (seps <= aperture)
        outcat = catalog.copy(deep=True)
        print(outcat)
        outcat.drop(outcat[~in_aperture_mask].index, inplace=True)
        print(outcat)
        obj_point_coords = list(zip(outcat.ra, outcat.dec))
        obj_points = [geometry.Point(point_coord) for point_coord in obj_point_coords]
        covered_mask = np.array([False]*len(outcat))
        for _, mask_object in regmask.iterrows():
            if mask_object['shape'] == 'circle':
                circle_center = geometry.point.Point(mask_object.ra, mask_object.dec)
                object_polygon = circle_center.buffer(mask_object.ax1.value)
            elif mask_object['shape'] == 'box':
                center_ra, center_dec = mask_object.ra, mask_object.dec
                dx, dy = mask_object.ax1.value, mask_object.ax2.value
                x = (center_ra - dx/2, center_ra - dx/2, center_ra + dx/2, center_ra + dx/2)
                y = (center_dec - dy/2, center_dec + dy/2, center_dec + dy/2, center_dec - dy/2)
                points = list(zip(x, y))
                object_polygon = geometry.Polygon(points)
            within = [point.within(object_polygon) for point in obj_points]
            covered_mask = covered_mask | np.array(within)
        
        return outcat[covered_mask]

    def tile(self, aperture_size = 120 * u.arcmin):
        #Tiles a region with the minimum number of aperture_size radius
        #circular sub-regions that covers the whole region
        #Note, this assumes a rectangular region
        bounds = self._region.bounds 
        size = aperture_size.to(u.degree).value
        dx = np.abs(bounds[0] - bounds[2])
        dy = np.abs(bounds[1] - bounds[3])
        num_x = math.ceil(dx/size) + 1
        num_y = math.ceil(dy/size) + 1
        bot_left_x, bot_left_y = min(bounds[0], bounds[2]), min(bounds[1], bounds[3])
        x_coords = [bot_left_x + i*size for i in range(num_x)]
        y_coords = [bot_left_y + i*size for i in range(num_y)]
        self._x_tiles = x_coords
        self._y_tiles = y_coords
        self._aperture_size = aperture_size
        self._tiled = True

    def load_bright_masks(self, path, type = 'HSC'):
        self._mask_path = path
        self._mask_type = type
        if self._tiled:
            if type == 'HSC':
                self._apply_hsc_mask(path)

    def get_bright_mask(self, center, aperture = 120*u.arcsec):
        start = time.time()

        regions = self.get_inside(center, 2*aperture)
        frames = []
        for region, data in regions.items():
            for id, subregion in data['subregions'].items():
                frames.append(subregion['ref'].get_mask_objects())
        
        full_frame = pd.concat(frames)
        coords = SkyCoord(full_frame['ra'], full_frame['dec'], unit="deg")
        distances = center.separation(coords)
        mask = (distances <= 2*aperture)
        return(full_frame[mask])

    def _apply_hsc_mask(self, path, filter = 'i'):
        contents = os.listdir(path)
        sel_reg_files = {}

        for item in contents:
            try: reg_int = int(item)
            except: continue
            try:
                current_region = self._subregions[reg_int]
                sel_reg_files.update({reg_int: {'patches': {}}})
            except:
                continue
            all_reg_files = [file for file in os.listdir(os.path.join(path, item))]
            for file in all_reg_files:
                data = file.split('-')
                if data[-1] == '{}.reg'.format(filter.upper()):
                    patch = data[2]
                    patch_key = tuple(map(int, patch.split(',')))
                    file_path = os.path.join(path, item, file)
                    data = self._parse_hsc_region_mask(file_path)
                    sel_reg_files[reg_int]['patches'].update({patch_key: data})
            current_region.add_mask_data(sel_reg_files[reg_int]['patches'])

    def _parse_hsc_region_mask(self, path):
        regions_list = []
        labels = ['shape', 'ra', 'dec', 'ax1', 'ax2', 'angle']

        with open(path) as f:
            for line in f:
                if line.startswith('circle'):
                    series = self._parse_circle_region(line)
                    regions_list.append(series)
                if line.startswith('box'):
                    series = self._parse_box_region(line)
                    regions_list.append(series)
        return pd.DataFrame(regions_list, columns = labels)

    def _parse_circle_region(self, line):
        shape = line.split("#")[0].strip()
        nums = [float(num) for num in re.findall(r'[+-]?\d+\.*\d*', line)]
        center = (nums[0], nums[1])
        radius = nums[2]*u.degree
        return ["circle", *center, radius, radius, 0]



    def _parse_box_region(self, line):
        shape = line.split("#")[0].strip()
        nums = [float(num) for num in re.findall(r'[+-]?\d+\.*\d*', line)]
        center = (nums[0], nums[1])
        r1, r2 = nums[2]*u.degree, nums[3]*u.deg
        angle = nums[4]
        if angle != 0:
            print("ANGLE: {}".format(angle))
        return ["box", *center, r1, r2, angle]

class region_old:
    def __init__(self, data):
        try: 
            self._type = data['type']
        except:
            self._type = 'None'

        self._center = SkyCoord(*data['center'], unit="deg")
        corner_ra, corner_dec = list(map(list, zip(*data['corners'])))
        self._corners = SkyCoord(corner_ra, corner_dec, unit="deg")
        if 'subregions' in data.keys():
            self._subregions = {key: region(value) for key, value in data['subregions'].items()}

    def get_inside(self, aperture_poly):
        subregions = {}
        if not hasattr(self, "_poly"):
            self._poly = geometry.Polygon(list(zip(self._corners.ra.degree, self._corners.dec.degree)))
        intersects = aperture_poly.intersects(self._poly)
        if not intersects:
            return False
        
        elif hasattr(self, "_subregions"):
            for id, subregion in self._subregions.items():
                inside = subregion.get_inside(aperture_poly)
                if inside or inside is None:
                    subregions.update({id: inside})
        return {'ref': self, 'subregions': subregions}
    
    def set_data(self, data):
        self._data = data
    def add_mask_data(self, maskdata):
        if type(maskdata) is pd.DataFrame:
            self._maskdata = maskdata
        elif hasattr(self, '_subregions'):
            for key, data in maskdata.items():
                try:
                    self._subregions[key].add_mask_data(data)
                except:
                    print("ERROR: Unable to add mask data for subregion {}".format(key))
    

    def get_mask_objects(self):
        if hasattr(self, '_maskdata'):
            return self._maskdata


    def attach_data(self):
        pass


            
if __name__ == '__main__':

    tracts = HSCSurvey.read_hsc_tracts('tracts_patches_W-w02.txt')
