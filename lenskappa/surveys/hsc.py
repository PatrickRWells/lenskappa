from numpy.lib.polynomial import poly
import re
import shapely
from lenskappa import region
from lenskappa import starmask
from lenskappa.surveys.survey import Survey, SurveyDataManager
from lenskappa.region import CircularSkyRegion, SkyRegion
from lenskappa.catalog import SkyCatalog2D, require_points, require_validation
from lenskappa.starmask import StarMaskCollection, RegStarMask, StarMask
import threading
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
from abc import ABCMeta, abstractmethod
from concurrent import futures
import traceback
import astropy.units as u



class HSCSurvey(Survey):
    
    datamanager = SurveyDataManager("hsc")
    def __init__(self, field, *args, **kwargs):
        self._field = field
        self._exec = futures.ProcessPoolExecutor()
        super().__init__(*args, **kwargs)
        #setup is called automatically by the superclass constructor

    def setup(self, *args, **kwargs):
        self._load_tract_data()
        self._load_starmasks(*args, **kwargs)

        self._load_catalog_data(*args, **kwargs)
        self._init_tracts(*args, **kwargs)

    def get_objects(self, region, masked=True, *args, **kwargs):
        """
        Get objects in a particular region.
        Parameters:
            region: shapely object defining the region of interest
            masked: whether to apply the bright star masks to the catalog first.
        
        """
        tracts = []
        patches = {}
        for id, tract in self._tracts.items():
            if tract.intersects(region):
                tracts.append(id)
                overlap_patches = tract.get_subregion_intersections(region)
                patches.update({id: overlap_patches})
        if not tracts:
            logging.error("Region deos not overlap with any tracts.")

        catalog = self._catalog.filter_by_subregions(patches)
        unmasked_catalog = catalog.get_objects_in_region(region)

        if masked:
            return self._starmasks.mask_catalog(unmasked_catalog, region, patches)
        else:
            return unmasked_catalog

    def mask_external_catalog(*args, **kwargs):
        pass

    
    def _load_tract_data(self, *args, **kwargs):
        """
        Loads data about the tracts and patches for the given HSC field
        Used to split up data, making it easier to manage
        
        """
        tractfile = self.datamanager.get_support_file_location({'type': 'tracts_patches', 'id': self._field})
        self._tract_masks = {}
        try:
            tracts = HSCSurvey._parse_tractfile(tractfile)
            self._tracts = HSCSurvey._parse_tractdata(tracts)
        except Exception as e:
            logging.error("Unable to read in tract data for HSC field {}".format(self._field))
            traceback.print_tb(e.__traceback__)
            return

    def _get_tract_mask(self, id):
        try:
            mask = self._tract_masks[id]
        except:
            mask = self._catalog['tract'] == id
            self._tract_masks.update({id: mask})
        return mask

    def _init_tracts(self, *args, **kwargs):
        for tract, region in self._tracts.items():
            self._catalog.add_subregion(tract, region)

    def _load_catalog_data(self, *args, **kwargs):
        """
        Returns catalog for a given HSC field
        """
        file = self.datamanager.get_file_location({'field': self._field, 'datatype': 'catalog'})
        self._cached_catalogs = {}
        if file:
            print("Reading in catalog for HSC field {}".format(self._field))

            catalog = pd.read_csv(file)
            self._catalog = hsc_catalog(catalog)
        else:
            logging.error("Unable to locate catalog file for HSC field {}".format(self._field))

    def _load_starmasks(self, *args, **kwargs):
        """
        Loads the bright star masks for the given field.
        Doesn't actually read the files in until they are needed
        Because the astropy Regions package is quite slow
        """
        logging.info("Reading in star masks for HSC field {}".format(self._field))
        mask_location = self.datamanager.get_file_location({'field': 'global', 'datatype': 'starmask'})
        paths = {}
        for tract in self._tracts.keys():
            try:
                wanted = tract in kwargs['tracts']
                if not wanted:
                    continue
            except:
                continue
            mask_path = os.path.join(mask_location, str(tract))
            if not os.path.exists(mask_path):
                logging.warnning("Unable to find bright star masks for HSC tract {}".format(tract))
            else:
                paths.update({tract: mask_path})
        try:
            band = kwargs['band']
        except:
            band = 'I'
        starmasks = {}
        for tract, path in paths.items():
            patch_files = [os.path.join(path, file) for file in os.listdir(path) if file.endswith(band + '.reg')]
            f = self._exec.submit(self._read_patch_maskfiles, patch_files)
            starmasks.update({tract: f})
        
        self._starmasks = hsc_mask(starmasks)


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
    def _get_patch_from_filename(filename):
        fname = os.path.basename(filename)
        patch = re.search(r"[0-9],[0-9]", fname)
        patch = patch.group()
        data = tuple(int(p) for p in patch.split(','))
        if len(data) == 2:
            return data        

    @staticmethod
    def _patch_tuple_to_int(patch_tuple):
        """
        Takes in a patch ID as a tuple and turns it into an int.
        This int can be used to look up objects in the catalof
        """
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

    
    @staticmethod
    def _parse_field_mask(maskdata):
        pass

    @staticmethod     
    def _read_patch_maskfiles(paths):

        patchdata = {}
        for file in paths:
            try:
                starmask_obj = RegStarMask.from_file(file, None)
                patch_tuple = HSCSurvey._get_patch_from_filename(file)
                patch_name = HSCSurvey._patch_tuple_to_int(patch_tuple)
                patchdata.update({patch_name: starmask_obj})
            except Exception as e:
                print(e)
                logging.warning("Couldn't find starmask at {}".format(file))
                
        return patchdata

    

class hsc_catalog(SkyCatalog2D):
    
    def __init__(self, cat):
        try:
            self._tract_names = cat['tract'].unique()
        except:
            logging.warning("Unable to initialize catalog tracts")
        super().__init__(cat)

    def _init_subregion(self, name, *args, **kwargs):
        if name in self._cat['tract'].unique():
            mask = self._cat['tract'] == name
            try:
                self._subregions[name]['mask'] = mask
                return mask
            except:
                logging.error("Unable to initialize mask for subregion {}".format(name))
        else:
            return super()._init_subregion(name)
    def filter_by_subregions(self, subregions):
        if type(subregions) == list:
            return super().filter_by_subregions(subregions)
        elif type(subregions) == dict:
            conditions = []
            for tract_id, patches in subregions.items():
                condition = {'tract': [tract_id], 'patch': patches}
                conditions.append(condition)
            return super().filter_by_columns(conditions)

class hsc_mask(StarMaskCollection):


    def __init__(self, masks):
        super().__init__(masks)
    
    

    def mask_catalog(self, catalog, region, patches):
        all_masks = []
        for tract_id, patches in patches.items():
            try:
                tract_mask = self._masks[tract_id]
            except Exception as e:
                logging.warning("No masks found for tract {}".format(tract_id))
                return
            
            try:
                mask = tract_mask.result()
                self._masks[tract_id] = mask
            except:
                mask = tract_mask
            
            for patch_id in patches:
                try:
                    all_masks.append(mask[patch_id])
                except Exception as e:
                    print(e)
                    logging.warning("Unable to find bright star mask for tract {} patch {}".format(tract_id, patch_id))

        final_mask = np.array([True]*len(catalog))

        for mask in all_masks:

            is_masked = mask.get_bool_region_mask(catalog, region)
            final_mask = final_mask & is_masked

        return catalog.from_dataframe(catalog[final_mask])            
            
    
    def mask_external_catalog(self, catalog, region):
        pass
    
