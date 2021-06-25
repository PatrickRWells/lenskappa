import parser
import numpy as np
from lenskappa.catalog import Catalog
from lenskappa.region import Region
from lenskappa.weighting import weighting
import logging
import astropy.units as u

class Counter:
    def __init__(self, field_catalog, control_catalog, control_region, aperture, control_mask = None, field_mask = None, *args, **kwargs):
        self._field_catalog = field_catalog
        self._control_catalog = control_catalog
        self._control_region = control_region
        self._aperture = aperture
        self._control_mask = control_mask
        self._validate_all()
        """
        Class for running the weighted number counts.
        Requires three inputs
        lens_catalog: <Catalog> A catalog for the field of interest
        field_catalog: <Catalog> A catalog for the control field
        field_region: <Region> A region defining the control field
        aperture: <float*astropy.u> Radius of the aperture to look at
        """
    def _validate_all(self):
        self._valid = True
        if not isinstance(self._field_catalog, Catalog):
            self._valid = False
            logging.error("Expected a catalog.Catalog object for the lens catalog")
        if not isinstance(self._control_catalog, Catalog):
            self._valid = False
            logging.error("Expected a catalog.Catalog object for the field catalog")
        if not isinstance(self._control_region, Region):
            self._valid = False
            logging.error("Expected a region.Region object for the field region")
        try:
            self._aperture = self._aperture.to(u.deg)
        except:
            self._valid = False
            logging.error("Expected an astropy Quantityt object for the aperture")

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
            results = self._get_weight_values(num_samples)
        else:
            results = self._delegate_weight_values(num_samples, threads)
        
        return results
        
    def _get_weight_values(self, num_samples, mutex = None, *args, **kwargs):
        
        for _ in range(num_samples):
            region = self._control_region.generate_tile(self._aperture)
            control_catalog = self._control_catalog.get_objects_in_region(region)
            field_catalog = self._field_catalog.get_objects_in_aperture(self._aperture)
            if self._control_mask is not None:
                control_catalog = self._control_mask.mask_catalog(control_catalog, region)
                field_catalog = self._control_mask.mask_external_catalog(field_catalog, region)
            
            field_weights = {key: weight.compute_weight(field_catalog) for key, weight in self._weights.items()}
            control_weights = {key: weight.compute_weight(control_catalog) for key, weight in self._weights.items()}

            row = self._parse_weight_vales(field_weights, control_weights)
            yield row
        

    def _delegate_weight_values(self, num_samples, threads):
        if threads <= 2:
            logging.warning("Minimum number of threads for a multithreaded run is 3"\
                            " (one supervisor thread, two worker threads.")
            return self._get_weight_values(num_samples)
        else:
            pass
        
    def _parse_weight_values(self, field_weights, control_weights):
        pass

if __name__ == '__main__':
    import toml
    config = toml.load('config/counts.toml')
    weight = count('zweight', config)