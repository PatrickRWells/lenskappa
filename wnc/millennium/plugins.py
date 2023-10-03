from dask.distributed.diagnostics.plugin import WorkerPlugin
from astropy.coordinates import SkyCoord
from heinlein import load_dataset
from cosmap.analysis.task import build_pipeline
from cosmap.analysis import utils
from cosmap.plugins import register

import astropy.units as u
from loguru import logger
import numpy as np
from itertools import product
from pydantic import BaseModel
import networkx as nx
import math
from functools import partial
from dask.distributed import get_worker
from heinlein import Region



@register
def generate_samples(sampler, *args, **kwargs):
    x_grid, y_grid = ms_sampler._generate_grid(sampler)
    #Assumption is that the aperture isn't changing during a run
    #maybe should change that
    coordinate_tuples = [ms_sampler.get_position_from_index(*index) for index in product(x_grid, y_grid)]
    coordinate_inputs = list(zip(*coordinate_tuples))
    ras = list(coordinate_inputs[0])
    decs = list(coordinate_inputs[1])
    coordinates = SkyCoord(ras, decs)
    return coordinates

@register
def initialize_sampler(sampler):
    pass


class ms_sampler:

    def _generate_grid(sampler):
        """
        Generates a grid of locations to compute weighted number counts on.
        Here, the locations are determined by the grid points defined in the
        millenium simulation.

        Params:

        aperture: Size of the aperture to consider. Should be an astropy quantity
        overlap: If 1, adjacent tiles do not overlap (centers have spacing of
            2*aperture). If above 1, tile spacing = 2*(aperture/overlap).
        """


        #First, find the corners of the tiling region.
        #Since we don't allow tiles to overlap with the edge of the field.

        spacing = sampler.analysis_parameters.ms_parameters.spacing

        aperture = max(sampler.sampler_parameters.sample_dimensions)
        min_pos = 0.0*u.degree + aperture
        max_pos = 4.0*u.degree - aperture
        bl_corner = ms_sampler.get_index_from_position(min_pos, min_pos)
        tr_corner = ms_sampler.get_index_from_position(max_pos, max_pos)
        #Since there's rounding involved above, check to make sure the tiles don't
        #Overlap with the edge of the field.
        min_pos_x, min_pos_y = ms_sampler.get_position_from_index(*bl_corner)
        max_pos_x, max_pos_y = ms_sampler.get_position_from_index(*tr_corner)

        
        min_vals = 0.0*u.degree
        max_vals = 4.0*u.degree
        pix_distance = 4.0*u.deg/4096.0

        x_diff = min_pos_x - min_vals
        y_diff = min_pos_y - min_vals
        x_index = bl_corner[0]
        y_index = bl_corner[1]

        #Make sure we're fully within the field
        if x_diff < aperture:
            x_index += 1
        if y_diff < aperture:
            y_index += 1
        bl_corner = (x_index, y_index)

        x_diff = max_vals - max_pos_x
        y_diff = max_vals - max_pos_y
        x_index = tr_corner[0]
        y_index = tr_corner[1]

        #Make sure we're fully within the field.
        if x_diff < aperture:
            x_index -= 1
        if y_diff < aperture:
            y_index -= 1
        tr_corner = (x_index, y_index)

        x_grid = np.arange(bl_corner[0], tr_corner[0], spacing)
        y_grid = np.arange(bl_corner[1], tr_corner[1], spacing)
        logger.info("With spacing of {}, millenium simulation grid has {} x {} = {} samples".format(spacing, len(x_grid), len(y_grid), len(x_grid)*len(y_grid)))        
        return x_grid, y_grid

    @staticmethod
    def get_position_from_index(x, y):
        """
        Returns an angular position (in radians) based on a given x, y index
        Where x,y are in the range [0, 4096]. This matches with the
        gridpoints defined by the millenium simulation.

        The position returned is with reference to the center of the field,
        so negative values are possible
        """
        l_field = (4.0*u.degree).to(u.radian) #each field is 4 deg x 4 deg
        n_pix = 4096.0
        l_pix = l_field/n_pix
        pos_x =  (x+0.5) * l_pix - 2.0*u.deg
        pos_y =  (y+0.5) * l_pix - 2.0*u.deg
        return pos_x,pos_y


    @staticmethod
    def get_index_from_position(pos_x, pos_y):
        """
        Returns the index of the nearest grid point given an angular position.
        
        """
        try:
            pos_x_rad = pos_x.to(u.radian)
            pos_y_rad = pos_y.to(u.radian)
        except:
            raise ValueError("Need angular distances to get kappa indices!")

        l_field = (4.0*u.degree).to(u.radian) #each field is 4 deg x 4 deg
        n_pix = 4096.0
        l_pix = l_field/n_pix

        x_pix = pos_x/l_pix - 0.5
        y_pix = pos_y/l_pix -0.5
        return int(round(x_pix.value)), int(round(y_pix.value))


@register
def generate_tasks(client, parameters: BaseModel, dependency_graph: nx.DiGraph, needed_dtypes: list, samples: list, chunk_size: int = 1000):
    fields_to_sample = getattr(parameters.analysis_parameters.ms_parameters, "fields_to_sample", None)

    pipeline_function = build_pipeline(parameters, dependency_graph)
    if fields_to_sample is not None:
        splits = [f.split("_") for f in fields_to_sample]
        field_tuples = [(int(s[0]), int(s[1])) for s in splits]
    else:
        field_numbers = list(range(8))
        field_tuples = list(product(field_numbers, field_numbers))

    n_workers = len(client.nthreads())
    sample_dimensions = parameters.sampling_parameters.sample_dimensions
    workers = list(client.has_what().keys())
    if n_workers > len(field_tuples):
        logger.warning("There are more workers than fields. This will result in some workers being idle.")
        n_workers = len(field_tuples)

    fields_per_worker = len(field_tuples) // n_workers
    extra_fields = len(field_tuples) % n_workers

    worker_fields = [field_tuples[i*fields_per_worker:(i+1)*fields_per_worker] for i in range(n_workers)]
    extra_fields = field_tuples[-extra_fields:]
    n_per_chunk = math.ceil(len(samples)/chunk_size)
    chunks = np.array_split(samples, n_per_chunk)

    for i in range(fields_per_worker):
        for chunk in chunks:
            tasks = []
            for j in range(n_workers):
                worker_field = worker_fields[j][i]
                f = partial(main_task, dtypes=needed_dtypes, sample_dimensions = max(sample_dimensions), pipeline_function = pipeline_function, field=worker_field)
                worker_tasks = client.submit(f, chunk, workers=[workers[j]], pure=False)
                tasks.append(worker_tasks)
            yield tasks

    for field in extra_fields:
        f = partial(main_task, dtypes=needed_dtypes, sample_dimensions = max(sample_dimensions), pipeline_function = pipeline_function, field=field)
        tasks = client.map(f, chunks, pure=False)
        yield tasks

   
def main_task(coordinates, sample_dimensions, dtypes, pipeline_function, field):
    worker = get_worker()
    dataset = worker.dataset
    field = tuple(int(f) for f in field)
    if (f := dataset.get_parameter("ms_field") is None) or f != field:
        dataset.set_field(field)
    results = []
    for i, c in enumerate(coordinates):
        if i% 100 == 0:
            logger.info(f"Worker {worker.address} has processed {i}/{len(coordinates)} samples")
        field_data = dataset.cone_search(c, sample_dimensions)
        sample_region = Region.circle(c, sample_dimensions)
        results.append(pipeline_function(field_data, sample_region))
    return results



    

def set_field(current_field, dask_worker):
    dataset = dask_worker.dataset
    dataset.set_field(current_field)
