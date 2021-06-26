from numpy.lib.polynomial import poly
import re
import shapely
from lenskappa import region
from lenskappa import starmask
from lenskappa.surveys.survey import Survey, SurveyDataManager
from lenskappa.region import SkyRegion
from lenskappa.catalog import SkyCatalog2D
from lenskappa.starmask import RegStarMask, StarMaskCollection, StarMask
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

class HSCSurvey(Survey):
    
    datamanager = SurveyDataManager("hsc")
    def __init__(self, field, *args, **kwargs):
        self._field = field
        super().__init__(*args, **kwargs)
        #setup is called automatically by the superclass constructor

    def setup(self, *args, **kwargs):
        self._load_tract_data()
        self._load_catalog_data(*args, **kwargs)
        self._load_starmasks(*args, **kwargs)        

    def get_objects(*args, **kwargs):
        pass

    def mask_external_catalog(*args, **kwargs):
        pass

    
    def _load_tract_data(self, *args, **kwargs):
        tractfile = self.datamanager.get_support_file_location({'type': 'tracts_patches', 'id': self._field})
        try:
            self._tracts = self.read_hsc_tracts(tractfile)
        except:
            logging.error("Unable to read in tract data for HSC field {}".format(self._field))
            return

    def _load_catalog_data(self, *args, **kwargs):
        """
        Returns catalog for a given HSC field
        """
        file = self.datamanager.get_file_location({'field': self._field, 'datatype': 'catalog'})
        if file:
            print("Reading in catalog for HSC field {}".format(self._field))
            data = pd.read_csv(file)
            cat = hsc_catalog(data)
            self._catalog =  cat
        else:
            logging.error("Unable to locate catalog file for HSC field {}".format(self._field))

    def _load_starmasks(self, *args, **kwargs):
        logging.info("Reading in star masks for HSC field {}".format(self._field))
        logging.info("This may take a while")

        mask_location = self.datamanager.get_file_location({'field': 'global', 'datatype': 'starmask'})
        paths = {}
        for tract in self._tracts.keys():
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
            starmasks.update({tract: lambda files=patch_files: self._read_patch_maskfiles(files)})
        self._starmasks = hsc_mask(starmasks)


    @staticmethod
    def read_hsc_tracts(tractfile):
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

    
    @staticmethod
    def _parse_field_mask(maskdata):
        pass
            
    def _read_patch_maskfiles(self, paths):

        patchdata = {}
        for file in paths:
            try:
                starmask_obj = RegStarMask.from_file(file, None)
                patch_name = self._get_patch_from_filename(file)
                patchdata.update({patch_name: starmask_obj})
            except Exception as e:
                print(e)
                logging.warning("Couldn't find starmask at {}".format(file))
                
        return patchdata


    def load_tracts(self, field):
        tractfile = self.datamanager.get_support_file_location({'type': 'tracts_patches', 'id': field})
        try:
            self._tracts = self.read_hsc_tracts(tractfile)
        except:
            logging.error("Unable to read in tract data for HSC field {}".format(field))
            return
    
    def _get_patch_from_filename(self, filename):
        fname = os.path.basename(filename)
        patch = re.search(r"[0-9],[0-9]", fname)
        patch = patch.group()
        data = tuple(int(p) for p in patch.split(','))
        if len(data) == 2:
            return data        


class hsc_catalog(SkyCatalog2D):
    def __init__(self, cat):
        super().__init__(cat)
    
    def init_tracts(self, tracts):
        self._tracts = tracts
        try:
            tracts = self._cat.tract.unique()
        except:
            logging.error("The HSC catalog does not contain data about "\
                          "which tract each object is in!")
        self._tract_masks = {}
        for tract_id in tracts:
            mask = self._cat.tract == tract_id
            self._tract_masks.update({tract_id: mask})

    def get_objects_in_region(self, region):
        """
        Overload of the base method for HSC data.
        This basically just makes the computation simpler by first filtering by tract
        And then actually applying the region.
        """
        overlap = []
        for tract_id, tract_region in self._tracts():
            if tract_region.overlaps(region):
                overlap.append(tract_id)
        
        if len(overlap) == 0:
            logging.error("Input region did not overlap with any tracts in this field")
            return
        
        else:
            mask = self._tract_masks[overlap[0]]
            for tract_id in overlap[1:]:
                mask = mask | self._tract_masks[tract_id]
        
        return self[mask]

class hsc_mask(StarMaskCollection):
    def __init__(self, masks):
        self._masks = masks
    
    @abstractmethod
    def mask_catalog(self, catalog, region):
        pass
    
    def mask_external_catalog(self, catalog, region):
        pass

if __name__ == "__main__":
    hsc = HSCSurvey("W02")
