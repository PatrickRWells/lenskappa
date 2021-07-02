from shapely import geometry
from concurrent import futures

import pandas as pd
import numpy as np
import math
import os
import logging
import operator
import traceback
import re

from lenskappa.surveys.survey import Survey
from lenskappa import SurveyDataManager
from lenskappa.region import SkyRegion
from lenskappa.catalog import SkyCatalog2D
from lenskappa.starmask import StarMaskCollection, RegStarMask


class HSCSurvey(Survey):

    datamanager = SurveyDataManager("hsc")
    def __init__(self, field, frame=None, *args, **kwargs):
        self._field = field
        self._exec = futures.ProcessPoolExecutor()
        self._frame = frame
        super().__init__(*args, **kwargs)
        #setup is called automatically by the superclass constructor

    def setup(self, *args, **kwargs):
        self._load_tract_data(*args, **kwargs)
        self._load_starmasks(*args, **kwargs)
        self._load_catalog_data(*args, **kwargs)
        self._init_tracts(*args, **kwargs)
        tract_objs = list(self._tracts.values())
        if self._frame is not None:
            self._region = self._frame
        else:
            self._region = SkyRegion.union(tract_objs)

        try:
            bind = kwargs['bind']
            if bind:
                logging.info("Binding sky region to catalog")
                self._bind()
        except:
            pass

    def _bind(self):
        """
        Shrinks the survey region to remove space that is
        not filled by the catalog.
        Should really only be used for testing purposes on
        small regions
        """
        points = self._catalog.get_points()
        multipoint = geometry.MultiPoint(points.to_numpy())
        bounds = multipoint.bounds
        box = geometry.box(*bounds)
        self._region = SkyRegion(box.centroid, box)



    def _get_patch_overlaps(self, region):
        tracts = []
        patches = {}
        for id, tract in self._tracts.items():
            if tract.intersects(region):
                tracts.append(id)
                overlap_patches = tract.get_subregion_intersections(region)
                patches.update({id: overlap_patches})
        if not tracts:
            logging.error("Region deos not overlap with any tracts.")

        return patches

    def generate_cirular_tile(self, radius, *args, **kwargs):
        pass


    def get_objects(self, region, masked=True, *args, **kwargs):
        """
        Get objects in a particular region.
        Parameters:
            region: shapely object defining the region of interest
            masked: whether to apply the bright star masks to the catalog first.

        """
        patches = self._get_patch_overlaps(region)
        catalog = self._catalog.filter_by_subregions(patches)
        unmasked_catalog = catalog.get_objects_in_region(region, *args, **kwargs)
        if len(unmasked_catalog) == 0:
            logging.error("Returned an empty dataframe")
            return unmasked_catalog
        unmasked_catalog = self.check_frame(region, unmasked_catalog)
        if masked:
            return self._starmasks.mask_catalog(unmasked_catalog, region, patches)
        else:
            return unmasked_catalog

    def mask_external_catalog(self, external_catalog, external_region, internal_region, *args, **kwargs):
        if not ( isinstance(external_region, SkyRegion) and isinstance(internal_region, SkyRegion) ):
            logging.error("Expected SkyRegion objects in mask_external_catalog")
            return None
        internal_area = internal_region.area
        external_area = external_region.area
        ratio = external_area/internal_area
        if ratio > 1.0:
            ratio = 1/ratio

        if ratio <= 0.99:
            logging.warning("Trying to put together two regions that are not the same size!"\
                            "I will only be able to center them on each other")

        patches = self._get_patch_overlaps(internal_region)
        new_catalog = self._starmasks.mask_external_catalog(external_catalog, external_region, internal_region, patches=patches, *args, **kwargs)
        new_catalog = self.check_frame(internal_region, new_catalog)
        return new_catalog


    def _load_tract_data(self, *args, **kwargs):
        """
        Loads data about the tracts and patches for the given HSC field
        Used to split up data, making it easier to manage

        """
        tractfile = self.datamanager.get_support_file_location({'type': 'tracts_patches', 'id': self._field})
        self._tract_masks = {}
        try:
            tracts = HSCSurvey._parse_tractfile(tractfile)
            self._tracts = self._parse_tractdata(tracts, *args, **kwargs)
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
            print("Reading in catalog for HSC field {}. This may take a while".format(self._field))

            catalog = pd.read_csv(file)
            self._catalog = hsc_catalog(catalog, *args, **kwargs)
            print("Done reading in catalog.")
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
            mask_path = os.path.join(mask_location, str(tract))
            if not os.path.exists(mask_path):
                logging.warning("Unable to find bright star masks for HSC tract {}".format(tract))
            else:
                paths.update({tract: mask_path})
        try:
            band = kwargs['band']
        except:
            band = 'I'
        starmasks = {}
        for tract, path in paths.items():
            patch_files = [os.path.join(path, file) for file in os.listdir(path) if file.endswith(band + '.reg')]
            f = self._exec.submit(self._read_patch_maskfiles, patch_files, tract)
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

    def _parse_tractdata(self, tractdata, *args, **kwargs):
        output = {}
        try:
            wanted_tracts = kwargs['tracts']
        except:
            wanted_tracts = []
        for name, tract in tractdata.items():
            if wanted_tracts and name not in wanted_tracts:
                    continue
            corners = tract['corners']
            center = tract['center']
            polygon = HSCSurvey._parse_polygon_corners(center, corners)
            if self._frame is not None:
                if  self._frame._polygon.disjoint(polygon):
                    logging.info("Skipping tract {} based on input frame.".format(name))
                    continue
                else:
                    logging.info("Including tract {}".format(name))

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
    def _read_patch_maskfiles(paths, tract = None):

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
        if tract is not None:
            print("Finished loading bright star masks for tract {}".format(tract))
        return patchdata



class hsc_catalog(SkyCatalog2D):

    def __init__(self, cat, *args, **kwargs):
        try:
            self._tract_names = cat['tract'].unique()
        except:
            logging.warning("Unable to initialize catalog tracts")
        super().__init__(cat, *args, **kwargs)

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



    def mask_catalog(self, catalog, region, patches, *args, **kwargs):
        all_masks = self._get_mask_objects_by_patch(patches, *args, **kwargs)

        final_mask = np.array([True]*len(catalog))

        for mask in all_masks:

            is_masked = mask.get_bool_region_mask(catalog, region)
            final_mask = final_mask & is_masked

        return catalog.apply_catalog_mask(~final_mask)

    def _get_mask_objects_by_patch(self, patches, *args, **kwargs):
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

        return all_masks


    def mask_external_catalog(self, external_catalog, external_region, internal_region, *args, **kwargs):
        try:
            patches = kwargs['patches']
        except:
            logging.error("Expected a list of patches associated with this region")
            return

        external_center = external_region.center
        internal_center = internal_region.center
        catalog = external_catalog.rotate(external_center, internal_center)
        masked_catalog = self.mask_catalog(catalog, internal_region, patches)
        return masked_catalog
