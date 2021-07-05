from re import M
import traceback
import numpy as np
from shapely import geometry
import multiprocessing as mp
import os
import pandas as pd
import math
import logging
import astropy.units as u

from lenskappa.catalog import Catalog
from lenskappa.region import Region, SkyRegion
from lenskappa.weighting import weighting
from lenskappa.surveys.survey import Survey


class Counter:

    def __init__(self, field_catalog, survey, region, mask=True, field_mask = None, *args, **kwargs):
        self._field_catalog = field_catalog
        self._reference_survey = survey
        self._field_region = region
        self._mask = mask
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

    def add_catalog_filter(self, filter, name, filter_type = 'absolute', which='both'):
        allowed_which = ['control', 'field', 'both']
        allowed_types = ['absolute', 'periodic']

        if which not in allowed_which:
            logging.error("Paramter which must be one of {}".format(allowed_which))
            return

        if filter_type not in allowed_types:
            logging.error("Allowed filter types are {}".format(allowed_types))
            return
        
        if filter_type == 'absolute':
            try:
                filters = self._absolute_filters
            except:
                self._absolute_filters = {}
                filters = self._absolute_filters

        elif filter_type == 'periodic':
            try:
                filters = self._periodic_filters
            except:
                self._periodic_filters = {}
                filters = self._periodic_filters
        
        filters.update({name: {'filter': filter, 'which': which}})
        
        

    def apply_absolute_filters(self, catalog, cattype, *args, **kwargs):
        """
        Absolute filters are filters that are applied to the catalog(s) at the beginning of the run.
        Any changes the filter makes to the catalog will stay with the catalog throughout the run.
        Examples include filtering out objects past a certain redshift
        """
        try:
            filters = self._absolute_filters
            for filter in filters.values():
                if filter['which'] == 'both' or filter['which'] == cattype:
                    catalog = self.apply_filter(catalog, filter)
            return catalog
                
        except:
            return catalog


    def apply_periodic_filters(self, catalog, cattype, *args, **kwargs):
        """
        Periodic filters are applied to the catalog(s) at each iteration after they have been modified by the survey,
        but before they have been passed to the weighting code.
        Examples include replacing the distances to all objects with r < 10" with 10", to avoid
        objects close to the center dominating.
        """

        try:
            filters = self._periodic_filters

            for filter in filters.values():
                if filter['which'] == 'both' or filter['which'] == cattype:
                    catalog = self.apply_filter(catalog, filter)
            return catalog
                
        except Exception as e:
            traceback.print_exc()
            return catalog

    def apply_filter(self, catalog, filter):
        """
        Theres no actual implementation difference between absolute and periodic filters.
        """
        output_catalog = filter['filter'](catalog)
        return output_catalog



    def get_weight_ratios(self, weights, num_samples = 100, output_file = "output.csv", threads = 1, *args, **kwargs):
        """
        get the weighted count ratios

        Paramters:
            weights: [<str>] list of names of weights, as defined in weighting.weightfns, use 'all' for all weights
            num_samples: Number of control apertures to generate
            threads: Number of threads to run 
        """
        self._output_fname = output_file
        self._validate_all(*args, **kwargs)
        if not self._valid:
            exit()

        self._field_catalog = self.apply_absolute_filters(self._field_catalog, 'field')
        self._field_catalog.get_distances(self._field_region.skycoord[0], unit=u.arcsec)
        self._reference_survey.handle_catalog_filter(self.apply_absolute_filters, cattype='control')
        
        if weights == 'all':
            self._weightfns = weighting.load_all_weights()
        elif type(weights) == list:
            self._weightfns = weighting.load_some_weights(weights)
        
        weight_data = pd.DataFrame(columns=list(self._weightfns.keys()))

        if self._weightfns is None:
            return


        if threads > 1:
            logging.warning("Multithreading has not been fully implemented. Falling back to a single thread")
            
        for index, row in enumerate(self._get_weight_values(num_samples)):
            weight_data = weight_data.append(row, ignore_index=True)
            if index % (num_samples/10) == 0:
                print("Completed {} out of {} samples".format(index, num_samples))
                self._write_output(weight_data)

        self._write_output(weight_data)
        
    def _get_weight_values(self, num_samples, mutex = None, *args, **kwargs):
        
        for _ in range(num_samples):

            tile = self._reference_survey.generate_circular_tile(self._radius)
            control_catalog = self._reference_survey.get_objects(tile, masked=self._mask, get_distance=True, dist_units = u.arcsec)
        
            if self._mask:

                field_catalog = self._reference_survey.mask_external_catalog(self._field_catalog, self._field_region, tile)

            else:

                field_catalog = self._field_catalog

            control_catalog = self.apply_periodic_filters(control_catalog, 'control')
            field_catalog = self.apply_periodic_filters(field_catalog, 'field')
            field_weights = {key: weight.compute_weight(field_catalog) for key, weight in self._weightfns.items()}
            control_weights = {key: weight.compute_weight(control_catalog) for key, weight in self._weightfns.items()}

            row = self._parse_weight_values(field_weights, control_weights)
            yield row

    def _weight_worker(self, num_samples, queue, index, *args, **kwargs):

        weight_data = pd.DataFrame(columns=list(self._weightfns.keys()))
        for index, row in self._get_weight_values(num_samples, *args, **kwargs):
            weight_data = weight_data.append(row, ignore_index=True)
            if index % (num_samples/10) == 0:
                logging.info("Completed {} out of {} samples".format(index, num_samples))
                logging.info("Writing output for first {} samples".format(index))
                queue.put(weight_data)
                weight_data = pd.DataFrame(columns=list(self._weightfns.keys()))

        self._running[index] = False
        if len(weight_data != 0):
            queue.put(weight_data)



    def _write_output(self, weights):
        weights.to_csv(self._output_fname)

    def _delegate_weight_values(self, num_samples, threads):
        if threads <= 2:
            logging.warning("Minimum number of threads for a multithreaded run is 3"\
                            " (one supervisor thread, two worker threads.")
            return self._get_weight_values(num_samples)
        else:
            num_cores = mp.cpu_count()
            if num_cores < threads:
                logging.warning("You requested more cores than this machine has available. "\
                                "I will reduce the number of threads")
                num_threads = num_cores
            else:
                num_threads = threads
            
            numperthread = math.ceil(num_samples/num_threads)
            self._queues = [mp.Queue() for _ in range(num_threads -1)]
            self._processes = [mp.Process(target=self._weight_worker, args=(numperthread,self._queues[i], i) ) for i in range(num_threads -1) ]
            self._running = [True]*(num_threads - 1)
            for process in self._processes:
                process.start()
            self._listen()
        
    def _listen(self):
        """
        Collects results from worker threads and periodically writes them to output
        """
        output_frame = pd.DataFrame(columns=list(self._weightfns.keys()))
        while any(self._running):
            dataframes = [q.get(block=True) for q in self._queues]
            subframe = pd.concat(dataframes)
            output_frame = output_frame.append(subframe, ignore_index=True)
            self._write_output(output_frame)
        final = []
        for queue in self._queues:
            try:
                df = queue.get(block=False)
                final.append(df)
            except:
                pass
        final_frame = output_frame.append(pd.concat(final), ignore_index=True)
        self._write_output(final_frame)
        for process in self._processes:
            process.join()

    def _parse_weight_values(self, field_weights, control_weights):
        np.seterr(invalid='raise')
        return_weights = {}
        for weight_name, weight_values in field_weights.items():
            try:
                field_weight = float(weight_values)
                control_weight = float(control_weights[weight_name])
                try:
                    ratio = field_weight/control_weight
                except Exception as e:
                    return_weights.update({weight_name: -1})
                    continue

                return_weights.update({weight_name: ratio})
                continue
            except:
                pass

            if type(weight_values) == pd.Series:
                field_weight = np.sum(weight_values)
                control_weight = np.sum(control_weights[weight_name])
                try:
                    ratio = field_weight/control_weight
                except Exception as e:
                    return_weights.update({weight_name: -1})
                    continue


                return_weights.update({weight_name: ratio})
            
        return pd.Series(return_weights)
            