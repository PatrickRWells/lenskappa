from abc import abstractmethod, ABC
from re import M
import numpy as np
from shapely import geometry
import multiprocess as mp
import os
import pandas as pd
import math
import logging
import astropy.units as u
import atexit

from lenskappa.datasets.dataset import DataSet
from lenskappa.catalog.catalog import Catalog
from lenskappa.spatial import Region
from lenskappa.weighting import weighting
from lenskappa.datasets.surveys.survey import Survey
from lenskappa.utils.multithreading import MultiThreadObject


class Counter(ABC):

    def __init__(self, comparison_dataset: DataSet, mask=True, *args, **kwargs):
        self._reference_survey = comparison_dataset
        self._mask = mask
        """
        Base class for running the weighted number counts.
        Requires two inputs
        comparison_dataset: The dataset being compared to
        mask: Whether or not to apply masks to this dataster
        field_mask
        """
    
    @abstractmethod
    def _validate_all(self, *args, **kwargs):
        """
        Validate the inputs. Will depend on what kind of weighting
        is being done
        """
        pass

    @abstractmethod
    def get_weights(self, *args, **kwargs):
        """
        Main function to be called by user
        Should validate the inputs.
        """

    @abstractmethod
    def _get_weight_values(self, *args, **kwargs):
        """
        Should be a generator that returns weight values        
        """
        pass

    @abstractmethod
    def _parse_weight_values(self, *args, **kwargs):
        """
        Should parse weight values into a pandas row
        """
        pass

    def add_catalog_filter(self, filter, name, filter_type = 'periodic', catalogs = 'all'):
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
        allowed_types =  ['periodic']

        if filter_type not in allowed_types:
            logging.error("Allowed filter types are {}".format(allowed_types))
            return

        elif filter_type == 'periodic':
            try:
                filters = self._periodic_filters
            except:
                self._periodic_filters = {}
                filters = self._periodic_filters

        filters.update({name: {'filter': filter, 'which': catalogs}})


    def apply_periodic_filters(self, catalog, catname, *args, **kwargs):
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
                if filter['which'] == 'all' or catname in filter['which']:
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


    def _delegate_weight_values(self, num_samples, threads, *args, **kwargs):
        """
        Delegates the weighting to multiple threads.
        Details of actual weighting are handled by individual
        subclasses 
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






class RatioCounter(Counter):

    def __init__(self, field_catalog, comparison_dataset, region: Region, mask=True, field_mask = None, *args, **kwargs):
        self._field_catalog = field_catalog
        self._field_mask = field_mask
        self._field_region = region
        super().__init__(comparison_dataset, mask, *args, **kwargs)

    def __getstate__(self):
        # Added because python 3.8 has some weird issues passing objects to 
        # Subprocesses that already have subprocesses in them
        state = self.__dict__.copy()
        if '_processes' in state.keys():
            state['_processes'] = None
            state['_queues'] = None
        return state

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



    def get_weights(self, weights, num_samples = 100, output_file = "output.csv", threads = 1, meds=False, *args, **kwargs):
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
        self._field_catalog.get_distances(self._field_region.skycoord[0], unit=u.arcsec)

        #Since the survey is not a Catalog object, it has a handler for filters.

        #Load the weights
        if weights == 'all':
            self._weightfns = weighting.load_all_weights()
        elif type(weights) == list:
            self._weightfns = weighting.load_some_weights(weights)

        #Initialize the dataframe for storage
        if meds:
            weight_names = list(self._weightfns.keys())
            self._weight_names = []
            for name in weight_names:
                self._weight_names.append(name)
                self._weight_names.append('_'.join([name, 'meds']))
        else:
            self._weight_names = list(self._weightfns.keys())
        weight_data = pd.DataFrame(columns=self._weight_names)

        #If no weights were loaded, terminate
        if self._weightfns is None:
            return

        sample_param = self._field_catalog.has_samples()
        if sample_param:
            self._has_catalog_samples = True
            self._generate_catalog_samples(sample_param)
        else:
            self._has_catalog_samples = False

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
        skipped_reference = 0
        skipped_field = 0
        while loop_i < num_samples:

            tile = self._reference_survey.generate_circular_tile(self._radius, *args, **kwargs)
            control_catalog = self._reference_survey.get_objects(tile, masked=self._mask, get_distance=True, dist_units = u.arcsec)
        


            if len(control_catalog) == 0:
                #Sometimes the returned catalog will be empty, in which case
                #we reject the sample
                skipped_reference += 1
                logging.warning("Found no objets for tile centered at {}".format(tile.skycoord[0]))
                logging.warning("In this thread, {} of {} samples have failed for this reason".format(skipped_reference, loop_i+skipped_reference+skipped_field+1))
                continue
            if self._has_catalog_samples:
                field_catalog = [self._prep_field_catalog(cat, internal_region=tile) for cat in self._sampled_catalogs]
            else:
                field_catalog = [self._prep_field_catalog(self._field_catalog, internal_region=tile)]
            
            if np.all([empty for empty in map(lambda cat: len(cat) == 0, field_catalog)]):
                skipped_field += 1
                logging.warning("Found no objects in field catalog after masking")
                logging.warning("In this thread, {} of {} samples have failed for this reason".format(skipped_field, loop_i+skipped_reference+skipped_field+1))
                continue
                
            control_catalog = self.apply_periodic_filters(control_catalog, 'control')
            control_weights={}
            field_weights={}
            for name in self._weight_names:
                if 'meds' not in name:
                    field_weights.update({name: [self._weightfns[name].compute_weight(c) for c in field_catalog]})
                    control_weights.update({name: self._weightfns[name].compute_weight(control_catalog)})
                    
                else:
                    wname = name.split('_')[0]
                    #I'd much rather fold these into a single funtion call
                    field_weights.update({name: [self._weightfns[wname].compute_weight(c, meds=True) for c in field_catalog]})
                    control_weights.update({name: self._weightfns[wname].compute_weight(control_catalog, meds=True)})

            row = self._parse_weight_values(field_weights, control_weights)
            loop_i += 1
            yield row
    
    def _prep_field_catalog(self, cat, *args, **kwargs):
            if self._mask:
                field_catalog = self._reference_survey.mask_external_catalog(cat, external_region=self._field_region, *args, **kwargs)

            else:
                field_catalog = cat

            field_catalog = self.apply_periodic_filters(field_catalog, 'field')
            return field_catalog



    def _generate_catalog_samples(self, sample_param, *args, **kwargs):
        #At present, we can only handle one sampled param at a time
        par = sample_param[0]
        field_catalogs = self._field_catalog.generate_catalogs_from_samples(par)
        self._sampled_catalogs = field_catalogs

    def _parse_weight_values(self, field_weights, control_weights):

        """
        Parse the values returned by the weighting code and compute the ratio
        """
        np.seterr(invalid='raise')

        return_weights = {weight_name: np.zeros(len(weight_values)) for weight_name, weight_values in field_weights.items()}

        for weight_name, weight_values in field_weights.items():
            try:
                control_weight = float(control_weights[weight_name])
            except:
                control_weight = np.sum(control_weights[weight_name])

            for index, value in enumerate(weight_values):

                try:
                    field_weight = float(value)
                    try:
                        ratio = field_weight/control_weight
                        return_weights[weight_name][index] = ratio
                    except Exception as e:
                        """
                        This should result from an overflow
                        """
                        return_weights[weight_name][index] = -1
                    continue
                except:
                    field_weight = np.sum(value)
                    try:
                        ratio = field_weight/control_weight
                        return_weights[weight_name][index] = ratio
                    except Exception as e:
                        return_weights[weight_name][index] = -1
                        continue

        return pd.DataFrame(return_weights)
    


class SingleCounter(Counter):

    def __init__(self, dataset: DataSet, mask = False, *args, **kwargs):
        """
        This counter does not compute ratios, it just gets the values of the weights
        for the control dataset passed.
        """
        super().__init__(dataset, mask, *args, **kwargs)
    
    def _validate_all(self, *args, **kwargs):
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
        if not isinstance(self._reference_survey, DataSet):
            logging.error("Reference surve is not a DataSet! Exiting...")
            exit()


    def get_weights(self, weights = 'all', aperture = 120*u.arcsec, sample_type='grid', output_file = "output.csv", threads = 1, output_positions = False, *args, **kwargs):
        self._output_fname = output_file
        self._validate_all(*args, **kwargs)
        if weights == 'all':
            weight_fns = weighting.load_all_weights()
        elif type(weights) is list:
            weight_fns = weighting.load_some_weights(weights)

        if sample_type == 'grid':
            print("Starting weighting...")

            self._get_weights_on_grid(aperture, weight_fns, output_file, output_positions = output_positions, *args, **kwargs)
        
    def _get_weights_on_grid(self,aperture, weights, output_file, output_positions = False, *args, **kwargs):
        columns = list(weights.keys())
        if output_positions:
            columns.insert(0, "ra")
            columns.insert(1, "dec")
        df = pd.DataFrame(columns=weights.keys())
        for index, reg in enumerate(self._reference_survey.get_ciruclar_tile(aperture, *args, **kwargs)):
            cat = self._reference_survey.get_objects_in_region(reg)
            cat = self.apply_periodic_filters(cat, "reference")
            row = {col: "" for col in columns}
            for name, weightfn in weights.items():
                weight_vals = weightfn.compute_weight(cat)
                row[name] = weight_vals
            if output_positions:
                pos = reg.skycoord[0]
                row['ra'] = pos.ra.deg
                row['dec'] = pos.dec.deg


            row = self._parse_weight_values(row)
            df = df.append(row, ignore_index=True)
            if index % 100 == 0:
                print("Completed {}".format(index))
                df.to_csv(output_file)



    def _get_weight_values(self, *args, **kwargs):
        pass

    def _parse_weight_values(self, row, *args, **kwargs):
        ret_val = {}
        for name, val in row.items():
            try:
                w_val = np.sum(val)
                ret_val.update({name: w_val})
            except:
                ret_val.update({name: val})
        
        return ret_val
 