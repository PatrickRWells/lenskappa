from abc import abstractmethod, ABC
import numpy as np
from shapely import geometry
import multiprocess as mp
import os
import pandas as pd
import math
import logging
import astropy.units as u
import atexit
from copy import copy

from lenskappa.weighting import weighting
from lenskappa.utils.multithreading import MultiThreadObject
from lenskappa.catalog import rotate

from heinlein.dataset import Dataset
from heinlein.dtypes.catalog import Catalog
from heinlein.region import BaseRegion, Region

class Counter(ABC):

    def __init__(self, comparison_dataset: Dataset, comparison_region: BaseRegion, mask=True, *args, **kwargs):
        self._reference_survey = comparison_dataset
        self._mask = mask
        self._comparison_region = comparison_region
        self._weight_params = {}
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

            except KeyError:
                logging.error("Error: filter {} does not have a 'which' parameter".format(name))

            except:
                logging.error("Unable to apply filter {}".format(name))
        
        return catalog


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
        #MultiThreadObject runs a check to make sure the number of threads is valid
        #So get the actual number from it after setting just to be safe.

        #The default numpy RNG is not thread safe, so we have to create several
        #This needs a more elegant solution
        splits = self._split_comparison_region(num_workers)
        splits = [Region.polygon(s) for s in splits]
        numperthread = math.ceil(num_samples/(num_workers))
        self._queues = [mp.Queue() for _ in range(num_workers)]
        self._processes = [mp.Process(target=weight_worker, args=(numperthread, splits[i], self._queues[i], i, self) ) for i in range(num_workers) ]
        for process in self._processes:
            atexit.register(process.terminate)

        self._running = [True]*(num_workers)
        for process in self._processes:
            process.start()
        print("Starting weighting...")

        self._listen(num_samples)

    def _split_comparison_region(self, threads):
        shapely_region = self._comparison_region.sky_geometry
        #NEED TO DO SOME SORT OF CHECKING HERE
        bounds = shapely_region.bounds
        dx = bounds[2] - bounds[0]
        dy = bounds[3] - bounds[1]
        if dx > dy:
            cuts = [bounds[0] + dx/(threads)*i for i in range(threads+1)]
            boxes = [geometry.box(b, bounds[1], cuts[i+1], bounds[3]) for i, b in enumerate(cuts[:-1])]
        else:
            cuts = [bounds[1] + dy/(threads)*i for i in range(threads+1)]
            boxes = [geometry.box(bounds[0], b, bounds[2], cuts[i+1]) for i, b in enumerate(cuts[:-1])]
        return boxes

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
                frames.append(output_frame)
                output_frame = pd.concat(frames)
                print("Completed {} out of {} samples".format(len(output_frame), num_samples))
                print("Writing output for first {} samples".format(index))

                self._write_output(output_frame)

        for i, process in enumerate(self._processes):
            print(f"Joined process {i}")
            process.join()


    def _weight_worker(self, num_samples, region, queue, thread_num, *args, **kwargs):
        """
        Worker function suitable for use in multithreaded runs.
        Expects a queue to communicate with the supervisor thread.
        """
        weight_data = pd.DataFrame(columns=list(self._weightfns.keys()))
        for index, row in enumerate(self._get_weight_values(num_samples, thread_num = thread_num, globals=globals, *args, **kwargs)):
            weight_data = weight_data.append(row, ignore_index=True)
            if index and (index % (int(num_samples/10)) == 0):

                print("Thread {} completed {} samples".format(thread_num, index))
                print("Sending to supervisor...")
                queue.put(weight_data)
                weight_data = pd.DataFrame(columns=list(self._weightfns.keys()))


        queue.put(weight_data)
        queue.put("done")



    def _write_output(self, weights):
        """
        Writes output.
        I/O should eventually be its own module
        """
        weights.to_csv(self._output_fname, index=False)


    def add_weight_params(self, params):
        self._weight_params.update(params)
        

class RatioCounter(Counter):

    def __init__(self, field_data, comparison_dataset, lens_region: BaseRegion, comparison_region: BaseRegion, mask=True, *args, **kwargs):
        self._field_catalog = field_data['catalog']
        if mask:
            self._field_mask = field_data['mask']
        else:
            self._field_mask = None
        self._field_region = lens_region
        super().__init__(comparison_dataset, comparison_region, mask, *args, **kwargs)

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
        if not isinstance(self._reference_survey, Dataset):
            logging.error("Expected a Survey object for the control field")
            self._valid = False

        if not isinstance(self._field_region, BaseRegion):
            logging.error("Expected a BaseRegion object for the lens field region")
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
        self._field_center = self._field_region.coordinate
        self._radius = self._field_region.radius
        self._field_catalog = self._field_catalog.get_data_from_region(self._field_region)



    def get_weights(self, weights, num_samples = 100, output_file = "output.csv", threads = 1, meds=False, notification_fraction = 0.1, *args, **kwargs):
        """
        get the weighted count ratios.

        Paramters:
            weights: [<str>] list of names of weights, as defined in weighting.weightfns, use 'all' for all weights
            num_samples: Number of control apertures to generate
            threads: Number of threads to run
        """
        self._notification_fraction = notification_fraction
        if threads > 1:
            MultiThreadObject.set_num_threads(threads)
        self._output_fname = output_file
        self._validate_all(*args, **kwargs)
        if not self._valid:
            exit()
        
        coords = self._field_catalog.coords
        distances = self._field_center.separation(coords).to(u.arcsec)
        self._field_catalog['r'] = distances

        #Since the survey is not a Catalog object, it has a handler for filters.

        #Load the weights
        if weights == 'all':
            self._weightfns = weighting.load_all_weights(self._weight_params)
        elif type(weights) == list:
            self._weightfns = weighting.load_some_weights(weights, self._weight_params)

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
        if self._field_mask is not None:
            self._field_catalog = self._field_catalog[self._field_mask]

        #Handle multithreaded runs
        if threads > 1:
            self._delegate_weight_values(num_samples, threads)

        #If only using one thread, just run the weighting
        else:
            print("Starting weighting...")
            for index, row in enumerate(self._get_weight_values(num_samples)):
                weight_data = pd.concat([weight_data, row], ignore_index=True)
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

            tile = self._comparison_region.generate_circular_tile(self._radius, *args, **kwargs)
            control_data = self._reference_survey.get_data_from_region(tile, ["catalog", "mask"])
            if control_data is None:
                skipped_reference += 1
                logging.warning("Found no objets for tile centered at {}".format(tile.coordinate))
                logging.warning("In this thread, {} of {} samples have failed for this reason".format(skipped_reference, loop_i+skipped_reference+skipped_field+1))
                continue

            control_catalog = control_data["catalog"]
            control_mask = control_data["mask"]
            distances = tile.coordinate.separation(control_catalog.coords).to(u.arcsec)
            control_catalog['r'] = distances

        


            if len(control_catalog) == 0:
                #Sometimes the returned catalog will be empty, in which case
                #we reject the sample
                skipped_reference += 1
                logging.warning("Found no objets for tile centered at {}".format(tile.coordinate))
                logging.warning("In this thread, {} of {} samples have failed for this reason".format(skipped_reference, loop_i+skipped_reference+skipped_field+1))
                continue
            if control_mask is not None:
                try:
                    control_catalog = control_catalog[control_mask]
                except ValueError:
                    print("Masking failed!")
                    print(control_catalog)
                    print(tile)
                    print(control_catalog.coords)
                    skipped_reference += 1
                    continue

                if len(control_catalog) == 0:
                    skipped_reference += 1
                    continue
                field_catalog = rotate(self._field_catalog, self._field_center, tile.coordinate)
                field_catalog = field_catalog[control_mask]
            else:
                field_catalog = self._field_catalog

            if self._field_mask is not None:
                control_catalog = rotate(control_catalog, tile.coordinate, self._field_center)
                control_catalog = control_catalog[self._field_mask]
            


            if len(field_catalog) == 0:
                skipped_field += 1
                logging.warning("Found no objects in field catalog after masking")
                logging.warning("In this thread, {} of {} samples have failed for this reason".format(skipped_field, loop_i+skipped_reference+skipped_field+1))
                continue
            control_catalog = self.apply_periodic_filters(control_catalog, 'control')
            field_catalog = self.apply_periodic_filters(field_catalog, 'field')
            control_weights={}
            field_weights={}
            for name in self._weight_names:
                if 'meds' not in name:
                    field_weights.update({name: self._weightfns[name].compute_weight(field_catalog)})
                    control_weights.update({name: self._weightfns[name].compute_weight(control_catalog)})
                    
                else:
                    wname = name.split('_')[0]
                    #I'd much rather fold these into a single funtion call
                    field_weights.update({name: self._weightfns[wname].compute_weight(field_catalog, meds=True)})
                    control_weights.update({name: self._weightfns[wname].compute_weight(control_catalog, meds=True)})

            row = self._parse_weight_values(field_weights, control_weights)
            loop_i += 1
            yield row
    
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

        return_weights = {}

        for weight_name, weight_values in field_weights.items():
            try:
                control_weight = float(control_weights[weight_name])
            except:
                control_weight = np.sum(control_weights[weight_name])

            try:
                field_weight = float(weight_values)
            except:
                field_weight = np.sum(weight_values)

            try:
                ratio = field_weight/control_weight
            except:
                ratio = -1
            return_weights[weight_name] = [ratio]
                
        return pd.DataFrame.from_dict(return_weights)
    

def weight_worker(num_samples, region, queue, thread_num, counter, *args, **kwargs):
    counter._comparison_region = region
    notification_fraction = counter._notification_fraction
    weight_data = pd.DataFrame(columns=list(counter._weightfns.keys()))
    for index, row in enumerate(counter._get_weight_values(num_samples, thread_num = thread_num, *args, **kwargs)):
        weight_data = pd.concat([weight_data,row], ignore_index=True)
        if index and (index % (int(num_samples/10)) == 0):
            print("Thread {} completed {} samples".format(thread_num, index))
            print("Sending to supervisor...")
            queue.put(weight_data)
            weight_data = pd.DataFrame(columns=list(counter._weightfns.keys()))
    queue.put(weight_data)
    queue.put("done")


class SingleCounter(Counter):

    def __init__(self, dataset: Dataset, mask = False, *args, **kwargs):
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
        if not isinstance(self._reference_survey, Dataset):
            logging.error("Reference surve is not a DataSet! Exiting...")
            exit()


    def get_weights(self, weights = 'all', meds=False, aperture = 120*u.arcsec, sample_type='grid', output_file = "output.csv", threads = 1, output_positions = False, *args, **kwargs):
        self._output_fname = output_file
        self._validate_all(*args, **kwargs)
        if weights == 'all':
            self._weightfns = weighting.load_all_weights()
        elif type(weights) is list:
            self._weightfns = weighting.load_some_weights(weights)

        if meds:
            weight_names = list(self._weightfns.keys())
            self._weight_names = []
            for name in weight_names:
                self._weight_names.append(name)
                self._weight_names.append('_'.join([name, 'meds']))
        else:
            self._weight_names = list(self._weightfns.keys())

        if sample_type == 'grid':
            print("Starting weighting...")

            self._get_weights_on_grid(aperture, output_file, output_positions = output_positions, *args, **kwargs)
        
    def _get_weights_on_grid(self,aperture, output_file, output_positions = False, *args, **kwargs):
        samples = self._reference_survey.has_samples()
        if samples and len(samples) != 1:
            samples = False
        columns = copy(self._weight_names)
        if output_positions:
            columns.insert(0, "ra")
            columns.insert(1, "dec")
        df = pd.DataFrame(columns=columns)
        for index, reg in enumerate(self._reference_survey.get_ciruclar_tile(aperture, *args, **kwargs)):
            catalog = self._reference_survey.get_objects_in_region(reg)
            if samples:
                cat = catalog.generate_catalogs_from_samples(samples[0], *args, **kwargs)
            else:
                cat = [catalog]
            
            cat = [self.apply_periodic_filters(c, "reference") for c in cat]
            weights = {}
            for name in self._weight_names:
                if 'meds' not in name:
                    weights.update({name: [self._weightfns[name].compute_weight(c) for c in cat]})
                    
                else:
                    wname = name.split('_')[0]
                    #I'd much rather fold these into a single funtion call
                    weights.update({name: [self._weightfns[wname].compute_weight(c, meds=True) for c in cat]})
                
            pos = reg.skycoord[0]


            row = self._parse_weight_values(weights, pos, output_positions)
            df = df.append(row, ignore_index=True)
            if index % 100 == 0:
                print("Completed {}".format(index))
                df.to_csv(output_file, index=False)



    def _get_weight_values(self, *args, **kwargs):
        pass

    def _parse_weight_values(self, weights, pos, output_positions, *args, **kwargs):
        return_weights = {weight_name: np.zeros(len(weight_values)) for weight_name, weight_values in weights.items()}

        for weight_name, weight_values in weights.items():
            for index, value in enumerate(weight_values):

                try:
                    final_weight = float(value)
                except:
                    final_weight = np.sum(value)

                return_weights[weight_name][index] = final_weight
        output = pd.DataFrame(return_weights)
        if output_positions:
            output['ra'] = pos.ra
            output['dec'] = pos.dec

        return output