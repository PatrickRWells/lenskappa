from abc import abstractmethod, ABC
from concurrent.futures import thread
import multiprocessing as mp
from operator import ge
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
import time

from lenskappa.weighting import weighting
from lenskappa.utils.multithreading import MultiThreadObject
from lenskappa.catalog import rotate
from lenskappa import output

from heinlein.dataset import Dataset
from heinlein.dtypes.catalog import Catalog
from heinlein.region import BaseRegion, Region




class Counter(ABC):

    def __init__(self, comparison_dataset: Dataset, comparison_region: BaseRegion, mask=True, *args, **kwargs):
        self._reference_survey = comparison_dataset
        self._mask = mask
        self._comparison_region = comparison_region
        self._weight_params = {}
        self.tnum = 0
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
        for s in splits:
            overlaps = self._reference_survey.get_region_overlaps(s)
            big_regions = [s_.name.split("_")[0] for s_ in overlaps]
            n_overlaps = len(set(big_regions))
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
        n_notify = int(num_samples*self._notification_fraction)
        last_notification = 0
        while np.any(self._running):
            for index, q in enumerate(self._queues):
                if not self._running[index]:
                    continue

                val = q.get()
                if val == "done":
                    self._running[index] = False
                else:
                    self._output.take_output(val)

            l = len(self._output)
            if l > (last_notification + n_notify):
                print("Completed {} out of {} samples".format(l, num_samples))
                print("Writing output for first {} samples".format(l))
                self._output.write_output(index=False)
                last_notification = l

        self._output.write_output(index=False)
        for i, process in enumerate(self._processes):
            print(f"Joined process {i}")
            process.join()


    def _weight_worker(self, num_samples, region, queue, thread_num, *args, **kwargs):
        """
        Worker function suitable for use in multithreaded runs.
        Expects a queue to communicate with the supervisor thread.
        """
        self.tnum = thread_num
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
        Writes weight_worker
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

        columns = ["ra", "dec"] + self._weight_names
        self._output = output.csvOutput(output_file, output.weightParser, columns)

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
            for index, row in enumerate(self._get_weight_values(num_samples, *args, **kwargs)):
                self._output.take_output(row)
                if index and index % (num_samples/10) == 0:
                    print("Completed {} out of {} samples".format(index, num_samples))
                    self._output.write_output(index=False)

            self._output.write_output(index=False)

    def _get_weight_values(self, num_samples, thread_num = 0, *args, **kwargs):
        """
        Generator that yields the weights. This is done in a complicated way,
        in order to conserve memory.
        """
        
        #Generate all the tiles we're going to use in advance
        samples = self._comparison_region.generate_circular_tiles(self._radius, num_samples)

        #Sort the samples by the region (or regions) that they fall into
        samples = self.partition_samples(samples)
        print(f"Thread {thread_num} finished initializing...")
        counts = [p.count("/") + 1 for p in samples.keys()]
        singles = []

        #We start with samples that fall onto multiple regions in the
        #Reference survey 
        for i, (regs, s) in enumerate(samples.items()):
            if counts[i] == 1:
                continue
            regs_to_get = regs.split("/")
            for s_ in self._get_samples(s, regs_to_get):
                yield(s_)
            for reg_ in regs_to_get:
                #Now, get the samples that fall into a single one of the regions
                if reg_ in singles:
                    continue
                try:
                    samples_in_reg = samples[reg_]
                except KeyError:
                    continue
                singles.append(reg_)
                for s_ in self._get_samples(samples_in_reg, [reg_]):
                    yield(s_)
            #Now, we dump the data from those regions
            self._reference_survey.dump_all()
    
        #Now we go through all the samples that fall in a single survey region 
        #that were NOT covered before
        for i, (reg, s) in enumerate(samples.items()):
            if counts[i] != 1 or reg in singles:
                continue
            else:
                for s_ in self._get_samples(s, [reg]):
                    yield s_
                self._reference_survey.dump_all()
    
    def _get_samples(self, samples, super_regions):
        num_samples = len(samples)
        rand = np.random.default_rng()
        rand_low = 0
        rand_high = len(super_regions)
        def generate_extra_tile():
            randint = rand.integers(rand_low, rand_high)
            reg = super_regions[randint]
            reg = self._reference_survey.get_region_by_name(reg)
            return reg.generate_circular_tile(self._radius)
        loop_i = 0
        while loop_i < num_samples:

            tile = samples[loop_i]
            print(tile.center)
            control_data = self._reference_survey.get_data_from_region(tile, ["catalog", "mask"])
            if control_data is None or len(control_data["catalog"]) == 0:
                samples[loop_i] = generate_extra_tile()
                continue

            control_catalog = control_data["catalog"]
            control_mask = control_data["mask"]
            distances = tile.coordinate.separation(control_catalog.coords).to(u.arcsec)
            control_catalog['r'] = distances
            print(f"Control catalog starts with {len(control_catalog)}")
            print(f"Field catalog starts with {len(self._field_catalog)}")
            if len(control_catalog) == 0:
                #Sometimes the returned catalog will be empty, in which case
                #we reject the sample
                samples[loop_i] = generate_extra_tile()
                continue
            if control_mask is None:
                print(f"Couldn't find a mask in control field centered at: {tile.coordinate} ")
                samples[loop_i] = generate_extra_tile()
                continue
            if control_mask is not None:
                try:
                    control_catalog = control_catalog[control_mask]
                    if len(control_catalog) == 0:
                        raise ValueError
                except ValueError:
                    samples[loop_i] = generate_extra_tile()
                    continue
            

                field_catalog = rotate(self._field_catalog, self._field_center, tile.coordinate)
                len_i = len(field_catalog)
                field_catalog = field_catalog[control_mask]
                len_f = len(field_catalog)
                print(f"Control mask removed {len_i - len_f} from field catalog.")

            else:
                field_catalog = self._field_catalog

            if self._field_mask is not None:
                len_i = len(control_catalog)
                control_catalog = rotate(control_catalog, tile.coordinate, self._field_center)
                control_catalog = control_catalog[self._field_mask]
                len_f = len(control_catalog)
                print(f"Field mask removed {len_i - len_f} from control catalog.")

            if len(control_catalog) == 0 or len(field_catalog) == 0:
                samples[loop_i] = generate_extra_tile()
                continue

            control_catalog = self.apply_periodic_filters(control_catalog, 'control')
            field_catalog = self.apply_periodic_filters(field_catalog, 'field')
            if len(field_catalog) == 0 or len(control_catalog) == 0:
                samples[loop_i] = generate_extra_tile()
                continue

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

            loop_i += 1
            row = {'center': tile.coordinate, 'field_weights': field_weights, 'control_weights': control_weights}
            yield row

    def partition_samples(self, samples, *args, **kwargs):
        overlaps = self._reference_survey.get_overlapping_region_names(samples)
        partitions = {}
        for i, sample in enumerate(samples):
            overlap = overlaps[i]
            if len(overlap) == 1:
                okey = overlap[0]
            else:
                overlap.sort()
                okey = "/".join(overlap)


            if okey in partitions.keys():
                partitions[okey].append(sample)
            else:
                partitions[okey] = [sample]
        return partitions

    def _generate_catalog_samples(self, sample_param, *args, **kwargs):
        #At present, we can only handle one sampled param at a time
        par = sample_param[0]
        field_catalogs = self._field_catalog.generate_catalogs_from_samples(par)
        self._sampled_catalogs = field_catalogs
    

def weight_worker(num_samples, region, queue, thread_num, counter, *args, **kwargs):
    counter._comparison_region = region
    n_notify = int(num_samples*counter._notification_fraction)
    start = time.time()
    for index, row in enumerate(counter._get_weight_values(num_samples, thread_num = thread_num, *args, **kwargs)):
        if index and (index % n_notify == 0):
            end = time.time()
            t = round(end - start)
            print(f"Thread {thread_num} completed {index} of {num_samples} samples!")
            print(f"The most recent {n_notify} samples took {t} seconds")
            print(f"On average, this is about {round(t/n_notify,1)} seconds/sample")
            start = time.time()
        queue.put(row)

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
        if weights == 'all':
            self._weightfns = weighting.load_all_weights(self._weight_params)
        elif type(weights) is list:
            self._weightfns = weighting.load_some_weights(weights, self._weight_params)

        if meds:
            weight_names = list(self._weightfns.keys())
            self._weight_names = []
            for name in weight_names:
                self._weight_names.append(name)
                self._weight_names.append('_'.join([name, 'meds']))
        else:
            self._weight_names = list(self._weightfns.keys())

        output_cols = ["ra", "dec"] + self._weight_names
        self._output = output.csvOutput(path=output_file, parser=output.singleWeightParser, columns=output_cols)
        if sample_type == 'grid':
            print("Starting weighting...")

            self._get_weights_on_grid(aperture, *args, **kwargs)
        
    def _get_weights_on_grid(self, aperture, overlap=1, *args, **kwargs):

        grid = self._reference_survey.generate_grid(aperture, overlap)


        for index, coord in enumerate(grid):
            catalog = self._reference_survey.cone_search(coord, aperture)["catalog"]
            cat_coords = catalog.coords
            distances = coord.separation(cat_coords).to(u.arcsec)
            catalog['r'] = distances

            catalog = self.apply_periodic_filters(catalog, "reference")
            weights = {}
            for name in self._weight_names:
                if 'meds' not in name:
                    weights.update({name: self._weightfns[name].compute_weight(catalog)})
                    
                else:
                    wname = name.split('_')[0]
                    #I'd much rather fold these into a single funtion call
                    weights.update({name: self._weightfns[wname].compute_weight(catalog, meds=True)})


            self._output.take_output({"center": coord, "weights": weights})
            if index % 100 == 0:
                print("Completed {}".format(index))
                self._output.write_output()



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