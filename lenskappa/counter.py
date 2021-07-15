from re import M
import numpy as np
from shapely import geometry
import multiprocessing as mp
import os
import pandas as pd
import math
import logging
import astropy.units as u
import atexit


from lenskappa.catalog import Catalog
from lenskappa.region import Region
from lenskappa.weighting import weighting
from lenskappa.surveys.survey import Survey
from lenskappa.utils.multithreading import MultiThreadObject


class Counter:

    def __init__(self, field_catalog, survey, region, mask=True, field_mask = None, *args, **kwargs):
        self._field_catalog = field_catalog
        self._reference_survey = survey
        self._field_region = region
        self._mask = mask
        """
        Class for running the weighted number counts.
        Requires three inputs
        field_catalog: <Catalog> A catalog for the field of interest
        survey: <Survey> A survey object represeting the control dataset
        field_region: <Region> A region defining the lens field
        mask: <bool> Whether or not to apply survey masks to lens field. Default True
        field_mask: The bright star mask for the lens field. Optional
        """

    def _validate_all(self, *args, **kwargs):
        """
        Validates the inputs. Is NOT run on construction.
        """
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
        """
        Add a filter to the cataog(s)

        Params:
            filter: <Filter> the filter to be applied
            name: <str> Name for the filter
            filter_type: one of ["absolute", "periodic"]
                Absolute filters are applied before any weighting
                Periodic filters are applied at each weighting step
            which: one of ["field", "control", "both"]
                Which catalog the filter should be applied to

        """
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
        Examples include filtering out objects past a certain redshift.
        """
        try:
            filters = self._absolute_filters

        except AttributeError:
            return catalog

        for name, filter in filters.items():
            try:
                if filter['which'] == 'both' or filter['which'] == cattype:
                    catalog = self.apply_filter(catalog, filter)
                    return catalog

            except KeyError:
                logging.error("Error: filter {} does not have a 'which' parameter".format(name))

            except:
                logging.error("Unable to apply filter {}".format(name))

    def apply_periodic_filters(self, catalog, cattype, *args, **kwargs):
        """
        Periodic filters are applied to the catalog(s) generated by each sample before
        they are passed to the weighting code.
        Examples include replacing the distances to all objects with r < 10" with 10", to avoid
        objects close to the center dominating.
        """

        try:
            filters = self._periodic_filters

        except AttributeError:
            return catalog


        for name, filter in filters.items():
            try:
                if filter['which'] == 'both' or filter['which'] == cattype:
                    catalog = self.apply_filter(catalog, filter)
                    return catalog

            except KeyError:
                logging.error("Error: filter {} does not have a 'which' parameter".format(name))

            except:
                logging.error("Unable to apply filter {}".format(name))


    def apply_filter(self, catalog, filter):
        """
        Applies a filter to the catalog.
        Theres no actual implementation difference between absolute and periodic filters at the moment.
        """
        output_catalog = filter['filter'](catalog)
        return output_catalog



    def get_weight_ratios(self, weights, num_samples = 100, output_file = "output.csv", threads = 1, *args, **kwargs):
        """
        get the weighted count ratios.

        Paramters:
            weights: [<str>] list of names of weights, as defined in weighting.weightfns, use 'all' for all weights
            num_samples: Number of control apertures to generate
            threads: Number of threads to run
        """
        if threads > 1:
            MultiThreadObject.set_num_threads(threads)
        self._output_fname = output_file
        self._validate_all(*args, **kwargs)
        if not self._valid:
            exit()

        self._field_catalog = self.apply_absolute_filters(self._field_catalog, 'field')
        self._field_catalog.get_distances(self._field_region.skycoord[0], unit=u.arcsec)

        #Since the survey is not a Catalog object, it has a handler for filters.
        self._reference_survey.handle_catalog_filter(self.apply_absolute_filters, cattype='control')

        #Load the weights
        if weights == 'all':
            self._weightfns = weighting.load_all_weights()
        elif type(weights) == list:
            self._weightfns = weighting.load_some_weights(weights)

        #Initialize the dataframe for storage
        weight_data = pd.DataFrame(columns=list(self._weightfns.keys()))

        #If no weights were loaded, terminate
        if self._weightfns is None:
            return

        #Handle multithreaded runs
        if threads > 1:
            self._delegate_weight_values(num_samples, threads)

        #If only using one thread, just run the weighting
        else:
            print("Starting weighting...")
            for index, row in enumerate(self._get_weight_values(num_samples)):
                weight_data = weight_data.append(row, ignore_index=True)
                if index and index % (num_samples/10) == 0:
                    print("Completed {} out of {} samples".format(index, num_samples))
                    self._write_output(weight_data)

            self._write_output(weight_data)

    def _get_weight_values(self, num_samples, *args, **kwargs):
        """
        Generator that yields the weights
        """
        loop_i = 0
        skipped = 0
        while loop_i < num_samples:

            tile = self._reference_survey.generate_circular_tile(self._radius, *args, **kwargs)
            control_catalog = self._reference_survey.get_objects(tile, masked=self._mask, get_distance=True, dist_units = u.arcsec)

            if len(control_catalog) == 0:
                #Sometimes the returned catalog will be empty, in which case
                #we reject the sample
                skipped += 1
                logging.warning("Found no objets for tile centered at {}".format(tile.skycoord[0]))
                logging.warning("In this thread, {} of {} samples have failed for this reason".format(skipped, loop_i+skipped+1))
                continue

            if self._mask:
                #Mask
                field_catalog = self._reference_survey.mask_external_catalog(self._field_catalog, self._field_region, tile)

            else:

                field_catalog = self._field_catalog

            control_catalog = self.apply_periodic_filters(control_catalog, 'control')
            field_catalog = self.apply_periodic_filters(field_catalog, 'field')
            field_weights = {key: weight.compute_weight(field_catalog) for key, weight in self._weightfns.items()}
            control_weights = {key: weight.compute_weight(control_catalog) for key, weight in self._weightfns.items()}

            row = self._parse_weight_values(field_weights, control_weights)
            loop_i += 1
            yield row

    def _delegate_weight_values(self, num_samples, threads, *args, **kwargs):
        """
        Delegates the weighting to multiple threads.
        """

        if threads <= 2:
            #Eventually, I'd like to rewrite this so the main thread does some of the work
            logging.warning("Minimum number of threads for a multithreaded run is 3"\
                            " (one supervisor thread, two worker threads.")
            return self._get_weight_values(num_samples)

        MultiThreadObject.set_num_threads(threads)
        num_threads = MultiThreadObject.get_num_threads()
        num_workers = num_threads - 1
        self._reference_survey.wait_for_setup()
        #MultiThreadObject runs a check to make sure the number of threads is valid
        #So get the actual number from it after setting just to be safe.

        #The default numpy RNG is not thread safe, so we have to create several
        #This needs a more elegant solution

        numperthread = math.ceil(num_samples/(num_workers))
        self._queues = [mp.Queue() for _ in range(num_workers)]
        self._processes = [mp.Process(target=self._weight_worker, args=(numperthread, self._queues[i], i) ) for i in range(num_workers) ]

        for process in self._processes:
            atexit.register(process.terminate)

        self._running = [True]*(num_workers)
        for process in self._processes:
            process.start()
        print("Starting weighting...")

        self._listen(num_samples)

    def _listen(self, num_samples):
        """
        Collects results from worker threads and periodically writes them to output
        """
        output_frame = pd.DataFrame(columns=list(self._weightfns.keys()))
        while np.any(self._running):
            frames = []
            for index, q in enumerate(self._queues):
                if not self._running[index]:
                    continue

                val = q.get()
                if type(val) == pd.DataFrame:
                    frames.append(val)
                elif val == "done":
                    self._running[index] = False
            if frames:
                output_frame = output_frame.append(pd.concat(frames), ignore_index=True)
                logging.info("Completed {} out of {} samples".format(len(output_frame), num_samples))
                logging.info("Writing output for first {} samples".format(index))

                self._write_output(output_frame)

        for process in self._processes:
            process.join()


    def _weight_worker(self, num_samples, queue, thread_num, *args, **kwargs):
        """
        Worker function suitable for use in multithreaded runs.
        Expects a queue to communicate with the supervisor thread.
        """

        weight_data = pd.DataFrame(columns=list(self._weightfns.keys()))
        for index, row in enumerate(self._get_weight_values(num_samples, thread_num = thread_num, globals=globals, *args, **kwargs)):
            weight_data = weight_data.append(row, ignore_index=True)
            if index and (index % (int(num_samples/10)) == 0):

                logging.info("Thread {} completed {} samples".format(thread_num, index))
                logging.info("Sending to supervisor...")
                queue.put(weight_data)
                weight_data = pd.DataFrame(columns=list(self._weightfns.keys()))


        queue.put(weight_data)
        queue.put("done")



    def _write_output(self, weights):
        """
        Writes output.
        I/O should eventually be its own module
        """
        weights.to_csv(self._output_fname)




    def _parse_weight_values(self, field_weights, control_weights):

        """
        Parse the values returned by the weighting code and compute the ratio
        """
        np.seterr(invalid='raise')

        return_weights = {}

        for weight_name, weight_values in field_weights.items():
            try:
                field_weight = float(weight_values)
                control_weight = float(control_weights[weight_name])
                try:
                    ratio = field_weight/control_weight
                except Exception as e:
                    """
                    This should result from an overflow
                    """
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
