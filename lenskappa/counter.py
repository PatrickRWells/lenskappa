import parser
import numpy as np
from pandas.core.indexes import multi
from shapely import geometry
from shapely.geometry.base import exceptNull
from lenskappa.catalog import Catalog
from lenskappa.region import Region, SkyRegion
from lenskappa.weighting import weighting
from lenskappa.surveys.survey import Survey
import multiprocessing as mp
import os
import pandas as pd
import math

import logging
import astropy.units as u

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
                print("Writing output for first {} samples".format(index))
                self._write_output(weight_data)
                print("Done")

        self._write_output(weight_data)
        
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
            

if __name__ == '__main__':
    from lenskappa.surveys import hsc
    from lenskappa.region import CircularSkyRegion
    from lenskappa.catalog import SkyCatalog2D
    from astropy.coordinates import SkyCoord

    box = geometry.box(30, -6, 39, -2)
    region = SkyRegion(box.centroid, box)

    lens_parmap = {'m_gal': 'demp_sm', 'z_gal': 'demp_photoz_best', 'z_s': 1.523}
    field_parmap = {'m_gal': 'demp_sm', 'z_gal': 'demp_photoz_best', 'z_s': 1.523}
    lens_field = SkyCatalog2D.read_csv("/Users/patrick/Documents/Current/Research/LensEnv/0924/weighting/lens_cat.csv", parmap=lens_parmap)
    survey = hsc("W02", parmap=field_parmap, frame=region)
    lens_field.add_column_bound('i_cmodel_mag', max=24.)
    survey.add_column_bound(column='i_cmodel_mag', max=24.)
    lens_field.add_column_bound('z_gal', max=1.523)
    survey.add_column_bound('z_gal', max=1.523)
    

    aperture = 120*u.arcsec
    center = SkyCoord(141.23246, 2.32358, unit="deg")
    lens_region = CircularSkyRegion(center, aperture)
    
    counter = Counter(lens_field, survey, lens_region, mask=True)
    counter.get_weight_ratios("all", output_file="out.csv", num_samples=1000, threads= 4, overwrite=True)
