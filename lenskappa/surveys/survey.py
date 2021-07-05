from abc import ABCMeta, abstractmethod
from copy import copy


import os
import logging
import toml
import logging
import argparse
import hashlib
import numpy as np

import lenskappa
from lenskappa.region import SkyRegion


class Survey(metaclass=ABCMeta):
    """
    Individual surveys are essentially plugins. They must have
    certain attributes and methods, to ensure they behave correctly
    in the core code, but how they do those things is totally up to them

    A survey should consist of several things:
        A Region, defining where the survey is on the sky
        A Catalog, containing the objects
        Optionally, a bright star mask
        A SurveyDataManger, defined as a class attribute

        It's up to the individual class to decide when and how to load these.
        It must define a setup method for this

        It should also have several methods for use during the weighted number counting
            mask_external_catalog: Mask a catalog that is NOT part of the survey,
                based on some defined region WITHIN the survey
            get_objects: Get objects from the survey's catalog, based on a region inside the survey area.
                            Should apply the bright star mask if it exists.
    """
    def __init__(self, *args, **kwargs):
        if not hasattr(self, "datamanager"):
            logging.critical("No data manager found for the survey!")
            return None
        self.setup(*args, **kwargs)
        self._validate()

    def _validate(self):
        try:
            region = self._region
        except:
            logging.error("No region found for the survey")
        
        try:
            catalog = self._catalog
        except:
            logging.error("Now catalog found for this survey")
    
    def frame(self, region):
        """
        Sets a new region for the survey.
        Designed for cases where you want to restrict where the code looks
        For example, if part of your survey is not covered in all the bands you want
        """
        if not isinstance(region, SkyRegion):
            logging.error("Expected a sky region object for the frame!")
            return
        else:
            self._region = region

    def handle_catalog_filter(self, filter_fn, *args, **kwargs):
        """
        Passes filters through to the catalog
        
        Parameters:
        filter_fn: Fn that will apply the filter(s)
        """

        self._catalog = filter_fn(self._catalog, *args, **kwargs)


    @abstractmethod
    def setup(self, *args, **kwargs):
        pass

    def generate_circular_tile(self, radius):
        """
        This should probably be overridden for some kinds of 
        """

        return self._region.generate_circular_tile(radius)

    @abstractmethod
    def mask_external_catalog(self, external_catalog, external_catalog_region, internal_region, *args, **kwargs):
        """
        Apply the bright star mask for a region inside the survey to a catalog from outside the survey region.
        
        Parameters:
            external_catalog: <catalog.Catalog> The catalog for the external objects
            external_region: <region.SkyRegion> A region defining the location of the catalog catalog
            internal_region: <region.SkyRegion> A region defining the location inside the survey to get the masks from
        """
    
    @abstractmethod
    def get_objects(self, internal_region, mask = True, get_dist = True, *args, **kwargs):
        """
        Get objects within in a particular region of the survey.
        Either with or without masking objects near brigh stars.

        Parameters:
            internal_region <region.SkyRegion> Region inside the survey area to get objects for
            mask: <bool> Whether or not to mask out objects based on the bright star masks
            get_dist: <bool> Whether or not to add the distance from the center of the region into the catalog
        
        """
        pass

    def check_frame(self, region, catalog, *args, **kwargs):
        """
        Checks to see if any objects fall outside the defined survey region
        And removes them from the catalog if so.
        """
        if self._region.contains(region):
            
            #If the input region falls completely within the survey region
            #Just return the original caatalog.
            
            return catalog

        points = catalog.get_points()
        mask = np.array([self._region.contains(point) for point in points])
        newcat = catalog.apply_boolean_mask(mask)
        return newcat    