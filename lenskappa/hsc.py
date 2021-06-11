from astropy.coordinates.solar_system import get_moon
from erfa.core import aper
from numpy.core.numeric import full
from numpy.lib.function_base import cov
from shapely.geometry import polygon
from lenskappa.mask import mask
import re
import regions
from astropy.coordinates import SkyCoord
from astropy import wcs
import astropy.units as u
from shapely import geometry
import pandas as pd
import numpy as np
import math
import os
import sys
import time
import matplotlib.pyplot as plt
import multiprocessing



class hsc_catalog:
    def __init__(self, file, params ={}):
        self._data = pd.read_csv(file)
        for key in params.keys():
            filter_type = key.split('_')[0]
            filter_col = '_'.join([filter_type, 'col'])

            if 'max' in key:
                max_val = params['_'.join([filter_type, 'max'])]
                filter_mask = self._data[filter_col] > max_val
                self._data.drop(self._data[filter_mask].index, inplace=True)
                print("Removed items based where {} is above {}".format(filter_col, max_val))

            if 'min' in key:
                min_val = params['_'.join([filter_type, 'min'])]
                filter_mask = self._data[filter_col] < min_val
                self._data.drop(self._data[filter_mask].index, inplace=True)
                print("Removed items based where {} is below {}".format(filter_col, min_val))

        if 'z_filter' in params.keys():
            z_col = params['z_col']
            z_mask = self._data[z_col] > params['z_filter']
            self._data.drop(self._data[z_mask].index, inplace=True)
        print("Splitting data...")
        self.subsets = {id: self._data[self._data.tract == id] for id in self._data.tract.unique()}
        self._patch_data = {}

    
    def get_objects_within_patches(self, tracts):
        frames = []
        for id, data in tracts.items():
            patch_ints = [self._patch_tuple_to_int(key) for key in data['subregions'].keys()]
            mask = (self.subsets[id].patch).isin(patch_ints)
            frames.append(self.subsets[id][mask])
        frame = pd.concat(frames)

        frame = frame.groupby(frame.index).first()
        return frame

    def _patch_tuple_to_int(self, patch_tuple):
        if patch_tuple[0] == 0:
            return patch_tuple[1]
        else:
            return 100*patch_tuple[0] + patch_tuple[1]
    
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
        top_subregions = {}
        aperture_deg = aperture.to(u.deg).value
        circle_center = geometry.point.Point(center.ra.degree, center.dec.degree)
        aperture_poly = circle_center.buffer(aperture_deg)
        for id, region in self._subregions.items():
            inside = region.get_inside(aperture_poly)
            if inside or inside is None:
                top_subregions.update({id: inside})
        return top_subregions
    
    def get_all_weight_ratios(self,ext_catalog = None, ext_catalog_center = None, weights = None, hsc_params = {}, ext_params = {}, threads = 1):
      
        if threads == 1:
            return_weights = {weight: [] for weight in weights.keys()}
            #for i in range(len(self._x_tiles)):
            for i in range(len(self._x_tiles)):
                #for j in range(len(self._y_tiles)):
                for j in range(len(self._y_tiles)):
                    print("{}, {}".format(i, j))
                    weight_vals = self.get_weight_ratios(i, j, ext_catalog, ext_catalog_center, weights, hsc_params, ext_params)
                    for name, val in weight_vals.items():
                        return_weights[name].append(val)
            return return_weights
        else:
            return self._get_all_weight_ratios_delegate(ext_catalog, ext_catalog_center, weights, hsc_params, ext_params, threads)
    
    def _get_all_weight_ratios_delegate(self, ext_catalog = None, ext_catalog_center = None, weights = None, hsc_params = {}, ext_params = {}, threads = 2):
        print("Threading...")
        num_x = len(self._x_tiles)
        numperthread = math.floor(num_x/(threads - 1))
        processes = []
        q = multiprocessing.Queue()
        for i in range(threads - 1):
            if i == threads - 1:
                range_ = [numperthread*i, num_x]
            else:
                range_ = [numperthread*i, numperthread*(i+1)]

            process = multiprocessing.Process(target=self._get_weight_ratios_range, args=[range_, ext_catalog, ext_catalog_center, weights, hsc_params, ext_params, q])
            process.start()
            processes.append(process)
        for process in processes:
            process.join()
        
        return_values = [q.get() for process in processes]
        return pd.concat(return_values, ignore_index=True)

    def _get_weight_ratios_range(self, x_range = [], ext_catalog = None, ext_catalog_center = None, weights = None, hsc_params = {}, ext_params = {}, queue = None):
        return_weights = {weight: [] for weight in weights.keys()}

        if x_range:
            for x in range(x_range[0], x_range[1]):
                for y in range(len(self._y_tiles)):
                    print("{}, {}".format(x,y))
                    weight_vals = self.get_weight_ratios(x, y, ext_catalog, ext_catalog_center, weights, hsc_params, ext_params)
                    for name, val in weight_vals.items():
                        return_weights[name].append(val)
        
        outframe = pd.DataFrame(return_weights)
        if queue is not None:
            queue.put(outframe)
        else:
            return outframe
                        
                        



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
        
        if len(cat) == 0:
            return {weight_name: -1.0 for weight_name in weights.keys()}
        field_weights = {key: weight.compute_weight(ext_cat_masked, ext_params) for key, weight in weights.items()}
        control_weights = {key: weight.compute_weight(cat, hsc_params) for key, weight in weights.items()}
        weight_ratios = {key: 0.0 for key in weights.keys()}
        for key, val in field_weights.items():
            try: 
                val_float = float(val)
                weight_ratios[key] = val_float/float(control_weights[key])
            except:
                field_val = sum(val)
                control_val = sum(control_weights[key])
                weight_ratios[key] = field_val/control_val
        return weight_ratios

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
        outcat.drop(outcat[~in_aperture_mask].index, inplace=True)
        if len(outcat) == 0:
            return outcat
        #print("Removed {} of {} objects because they were too far away...".format(len(catalog) - len(outcat), len(catalog)))
        obj_point_coords = list(zip(outcat.ra, outcat.dec))
        obj_points = [geometry.Point(point_coord) for point_coord in obj_point_coords]
        covered_mask = np.array([False]*len(outcat))
        import matplotlib.pyplot as plt
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
        inside_top_region = [point.within(self._region) for point in obj_points]
        return outcat[(~covered_mask) & inside_top_region]

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
        print("Split region into {} tiles ({} by {})".format(num_x*num_y, num_x, num_y))

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

class region:
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


def read_hsc_tracts(tractfile):
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
            
if __name__ == '__main__':
    pass