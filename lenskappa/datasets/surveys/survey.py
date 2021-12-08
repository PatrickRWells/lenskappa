from abc import ABC, abstractmethod
from copy import copy


import os
import logging
import toml
import logging
import argparse
import hashlib
import numpy as np

import lenskappa
from lenskappa.spatial import SkyRegion
from lenskappa.datasets.dataset import SkyDataSet


class Survey(SkyDataSet):
    """
    A survey is a DataSet that contains information in a large region
    of the physical sky (i.e. not a simulation)

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
    def __init__(self, name, *args, **kwargs):

        super().__init__(name, *args, **kwargs)
        if not hasattr(self, "_datamanager"):
            logging.critical("No data manager found for the survey!")
            return None
        self.setup(*args, **kwargs)
        self._validate()

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
        Passes filters through to a catalog

        Parameters:
        filter_fn: Fn that will apply the filter(s)
        """

        self._catalog = filter_fn(self._catalog, *args, **kwargs)


    @abstractmethod
    def setup(self, *args, **kwargs):
        pass

    def generate_circular_tile(self, radius, *args, **kwargs):
        """
        This should probably be overridden for some kinds of
        """

        return self._region.generate_circular_tile(radius, *args, **kwargs)

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

    @abstractmethod
    def wait_for_setup(self, *args, **kwargs):
        """
        At present, surveys are only used to retrieve values during weighting.
        As such, there is no need for state to be shared be subprocesses when
        multiprocessing is added. This method will be called when running weighting in
        parallel, and should block execution until the survey object is ready to recieve
        requests.
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

