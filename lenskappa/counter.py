import parser
import numpy as np
from lenskappa.catalog import Catalog
from lenskappa.region import Region
from lenskappa.weighting import weighting
from lenskappa.surveys.survey import Survey
import os
import pandas as pd

import logging
import astropy.units as u
class Counter:

    def __init__(self, field_catalog, survey, region, mask=True, field_mask = None, output_file = "output.csv", *args, **kwargs):
        self._field_catalog = field_catalog
        self._reference_survey = survey
        self._field_region = region
        self._mask = mask
        self._output_fname = output_file
        self._validate_all(*args, **kwargs)
        """
        Class for running the weighted number counts.
        Requires three inputs
        lens_catalog: <Catalog> A catalog for the field of interest
        field_catalog: <Catalog> A catalog for the control field
        field_region: <Region> A region defining the control field
        aperture: <float*astropy.u> Radius of the aperture to look at
        """
    def _validate_all(self, *args, **kwargs):
        self._valid = True
        if not isinstance(self._field_catalog, Catalog):
            self._valid = False
            logging.error("Expected a catalog.Catalog object for the lens catalog")
        if not isinstance(self._reference_survey, Survey):
            logging.error("Expected a Survey object for the control field")
            self._valid = False

        if not isinstance(self._field_region, Region):
            logging.error("Expected a Region object for the lens field region")
            self._valid = False
        
        if os.path.exists(self._output_fname):
            logging.warning("File {} already exists.".format(self._output_fname))
            try:
                overwrite = kwargs['overwrite']
                if overwrite:
                    logging.warning("This file will be overwritten")
                else:
                    logging.warning("Set overwrite=True to overwrite this file anyway")
                    exit()
            except:
                logging.warning("Set overwrite=True to overwrite this file anyway")
                exit()

        


        #Remove all objects from the field catalog that fall outside the region
        self._field_center, self._radius = self._field_region.skycoord
        self._field_catalog = self._field_catalog.get_objects_in_region(self._field_region)


    def get_weight_ratios(self, weights, num_samples = 100, threads = 1, *args, **kwargs):
        """
        get the weighted count ratios

        Paramters:
            weights: [<str>] list of names of weights, as defined in weighting.weightfns, use 'all' for all weights
            num_samples: Number of control apertures to generate
            threads: Number of threads to run 
        """

        if weights == 'all':
            self._weightfns = weighting.load_all_weights()
        elif type(weights) == list:
            self._weightfns = weighting.load_some_weights(weights)

        if self._weightfns is None:
            return

        if threads == 1:
            
            for row in self._get_weight_values(num_samples):
                print(row)
        else:
            results = self._delegate_weight_values(num_samples, threads)
        
        return results
        
    def _get_weight_values(self, num_samples, mutex = None, *args, **kwargs):
        
        for _ in range(num_samples):

            tile = self._reference_survey.generate_circular_tile(self._radius)
            control_catalog = self._reference_survey.get_objects(tile, masked=self._mask, get_distance=True)
        
            if self._mask:

                field_catalog = self._reference_survey.mask_external_catalog(self._field_catalog, self._field_region, tile)

            else:

                field_catalog = self._field_catalog

            field_weights = {key: weight.compute_weight(field_catalog) for key, weight in self._weightfns.items()}
            control_weights = {key: weight.compute_weight(control_catalog) for key, weight in self._weightfns.items()}

            row = self._parse_weight_values(field_weights, control_weights)
            yield row
        

    def _delegate_weight_values(self, num_samples, threads):
        if threads <= 2:
            logging.warning("Minimum number of threads for a multithreaded run is 3"\
                            " (one supervisor thread, two worker threads.")
            return self._get_weight_values(num_samples)
        else:
            pass
        
    def _parse_weight_values(self, field_weights, control_weights):
        return_weights = {}
        for weight_name, weight_values in field_weights.items():
            try:
                field_weight = float(weight_values)
                control_weight = float(control_weights[weight_name])
                ratio = field_weight/control_weight
                return_weights.update({weight_name: ratio})
                continue
            except:
                pass

            if type(weight_values) == pd.Series:
                field_weight = np.sum(weight_values)
                control_weight = np.sum(control_weights[weight_name])
                ratio = field_weight/control_weight
                return_weights.update({weight_name: ratio})
            else:
                exit()
            
        return return_weights
            

if __name__ == '__main__':
    from lenskappa.surveys import hsc
    from lenskappa.region import CircularSkyRegion
    from lenskappa.catalog import SkyCatalog2D
    from astropy.coordinates import SkyCoord
    lens_parmap = {'m_gal': 'demp_sm', 'z_gal': 'demp_photoz_best', 'z_s': 1.523}
    field_parmap = {'m_gal': 'demp_sm', 'z_gal': 'demp_photoz_best', 'z_s': 1.523}
    lens_field = SkyCatalog2D.read_csv("/Users/patrick/Documents/Current/Research/LensEnv/0924/weighting/lens_cat.csv", parmap=lens_parmap)
    survey = hsc("W02_test", tracts=[9005, 9006], parmap=field_parmap)

    center = SkyCoord(141.23246, 2.32358, unit="deg")
    internal_center = SkyCoord(32, -2, unit="deg")
    aperture = 120*u.arcsec
    lens_region = CircularSkyRegion(center, aperture)
    counter = Counter(lens_field, survey, lens_region, mask=True, output_file="out.csv")
    counter.get_weight_ratios("all", num_samples=10)
