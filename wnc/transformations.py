from functools import partial
from loguru import logger
from typing import Callable
import astropy.units as u
from cosmap.analysis import CosmapBadSampleError
import json
from pathlib import Path
from dask.distributed import get_worker
from copy import copy
import numpy as np
import psutil

column_file = Path(__file__).parent / "columns.json"

with open(column_file, "r") as f:
    columns = json.load(f)


def build_filter(z_s, inner_radius, radius, limiting_magnitude, column_aliases) -> Callable:
    def filter_function(catalog, radius, limiting_magnitude):
        mask = catalog[column_aliases["mag_i"]] < limiting_magnitude
        mask & (catalog["radius"] > inner_radius)
        mask & (catalog["radius"] < radius)
        mask & (catalog[column_aliases["redshift"]] < z_s)
        return catalog[mask]
    return partial(filter_function, radius=radius, limiting_magnitude=limiting_magnitude)


class Setup:

    @staticmethod
    def build_filter_sets(lens_parameters, inner_radius, radius, limiting_magnitude, dataset_name) -> dict:
        if dataset_name not in columns.keys():
            raise ValueError(f"Dataset {dataset_name} not found in columns.json")
        column_aliases = columns[dataset_name]
        #make sure radius is a list
        if radius.isscalar:
            radius = [radius]
        #make sure limiting_magnitude is a list
        if isinstance(limiting_magnitude, float):
            limiting_magnitude = [limiting_magnitude]
        filters = {}
        #Create keys for the filter sets by combining the radius and limiting magnitude
        for lens_name, lens_data in lens_parameters.items():
            for r in radius:
                for m in limiting_magnitude:
                    #get limiting magnitude without any a period
                    mag_key = str(m).replace(".", "")
                    #drop any trailing zeros
                    mag_key = mag_key.rstrip("0")
                    #make a key by combing the radius and limiting magnitude
                    key = f"{lens_name}_{int(r.value)}_{mag_key}"
                    filter = build_filter(lens_data.source_redshift, inner_radius, r, m, column_aliases)
                    filters[key] = filter
        return filters
    
    @staticmethod
    def build_output_paths(filter_sets, base_output_path) -> dict:
        lenses = set([f.split("_")[0] for f in filter_sets.keys()])
        paths = {}
        for lens in lenses:
            if (base_output_path / lens).exists():
                logger.warning(f"Output directory {base_output_path / lens} already exists. Continuing anyway...")
            (base_output_path / lens).mkdir(parents=False, exist_ok=True)
            included_filters = [f for f in filter_sets.keys() if f.startswith(lens)]
            lens_paths = {f: base_output_path / lens / f"{f.replace(f'{lens}_', '')}.csv" for f in included_filters}
            paths.update(lens_paths)
        for path in paths.values():
            if path.exists():
                logger.warning(f"Output file {path} already exists. Continuing anyway...")
            path.parents[0].mkdir(parents=False, exist_ok=True)
        return paths

    @staticmethod
    def get_columns(weights):
        columns = ["ra", "dec"]
        for w in weights:
            columns.append(w)
            columns.append(f"{w}_meds")
        return columns  
      
class Main:

    @staticmethod
    def get_lens_data(radius, lens_parameters, dtypes, *args, **kwargs):
        try:
            rad = max(radius)
        except TypeError:
            rad = radius
        worker = get_worker()
        try:
            return worker.lens_data
        except AttributeError:
            logger.info("Loading lens data...")
            worker.lens_data = {}
            for lens_name, lens_data in lens_parameters.items():
                lens_coordinate = lens_data["lens_coordinate"]
                lens_data_ = worker.dataset.cone_search(lens_coordinate, rad, dtypes=dtypes)
                lens_distances = lens_coordinate.separation(lens_data_["catalog"]['coordinates'])
                lens_data_["catalog"]["radius"] = lens_distances.to(u.arcmin)
                worker.lens_data[lens_name] = lens_data_
                worker.dataset.dump_all()
            logger.info(f"Finished loading lens data")
            return worker.lens_data

    @staticmethod
    def apply_masks(catalog, mask, lens_data, lens_parameters, sample_region, *args, **kwargs):
        output = {}
        if len(catalog) == 0:
            raise CosmapBadSampleError
        distances = sample_region.coordinate.separation(catalog['coordinates'])
        catalog["radius"] = distances.to(u.arcmin)

        for lens_name, lens_data in lens_data.items():
            lens_coordinate = lens_parameters[lens_name]["lens_coordinate"]
            lens_catalog = lens_data["catalog"]
            lens_mask = lens_data["mask"]
            masked_lens_catalog = lens_mask.mask(lens_catalog)
            masked_field_catalog = mask.mask(catalog)

            if len(masked_lens_catalog) == 0 or len(masked_field_catalog) == 0:
                raise CosmapBadSampleError

            rotated_lens_catalog = rotate(masked_lens_catalog, lens_coordinate, sample_region.coordinate)
            rotated_field_catalog = rotate(masked_field_catalog, sample_region.coordinate, lens_coordinate)

            fully_masked_lens_catalog = mask.mask(rotated_lens_catalog)
            fully_masked_field_catalog = lens_mask.mask(rotated_field_catalog)

            if len(fully_masked_lens_catalog) == 0 or len(fully_masked_field_catalog) == 0:
                raise CosmapBadSampleError
            output.update({lens_name: {"lens_catalog": fully_masked_lens_catalog, "field_catalog": fully_masked_field_catalog}})

        return output

    @staticmethod
    def apply_filters(masked_catalogs, filters, output_paths, sample_region, *args, **kwargs):
        output = {} 
        for lens_name, catalogs in masked_catalogs.items():
            
            field_catalog = catalogs["field_catalog"]
            lens_catalog = catalogs["lens_catalog"]
            filters_for_lens = {name: filter for name, filter in filters.items() if name.startswith(lens_name)}
            field_catalogs = {name: filter(field_catalog) for name, filter in filters_for_lens.items()}
            lens_catalogs = {name: filter(lens_catalog) for name, filter in filters_for_lens.items()}
            if any([len(c) == 0 for c in field_catalogs.values()]):
                raise CosmapBadSampleError
            if any([len(c) == 0 for c in lens_catalogs.values()]):
                raise CosmapBadSampleError
            output.update({lens_name: {"lens_catalogs": lens_catalogs, "field_catalogs": field_catalogs}})
        return output


    @staticmethod
    def count(catalogs, sample_region, dataset_name, lens_parameters, *args, **kwargs):
        outputs = {}
        column_aliases = columns[dataset_name]
        redshift_alias = column_aliases["redshift"]
            

        for lens_name, lens_catalogs_ in catalogs.items():
            z_source = lens_parameters[lens_name]["source_redshift"]
            lens_catalogs = lens_catalogs_["lens_catalogs"]
            field_catalogs = lens_catalogs_["field_catalogs"]
            for name, catalog in lens_catalogs.items():
                lens_weights = compute_weightfns(catalog, redshift_alias, z_source)
                field_weights = compute_weightfns(field_catalogs[name], redshift_alias, z_source)


                final_weights = {name: lens_weights[name] / field_weights[name] for name in lens_weights.keys()}
                outputs[name] = {
                    "ra": sample_region.coordinate.ra.deg,
                    "dec": sample_region.coordinate.dec.deg,
                    **final_weights
                }

        return outputs



def compute_weightfns(catalog, redshift_key, source_redshift):
    weightfns = get_weightfns()
    lcat = len(catalog)
    weight_vals = {key: weightfn(catalog, redshift_key, source_redshift) for key, weightfn in weightfns.items()}
    weight_vals["z/r"] = weight_vals["z"] * weight_vals["1/r"]
    meds = {}
    for name, vals in weight_vals.items():
        meds[f"{name}_meds"] = np.median(vals)*lcat
    weight_vals = {n: np.sum(vals) for n, vals in weight_vals.items()}
    weight_vals.update(meds)
    return weight_vals

def get_weightfns():
    return {"n": ncount, "1/r": oneoverr, "z": zweight}

def rotate(catalog, center, new_center):
    original = catalog['coordinates']
    separations = original.separation(center)
    pas = original.position_angle(center)
    new_coords = new_center.directional_offset_by(pas - 180*u.deg, separations)
    new_cat = copy(catalog)
    new_cat['coordinates'] = new_coords
    return new_cat

def ncount(catalog, *args, **kwargs):
    return np.ones(len(catalog))

def oneoverr(catalog, *args, **kwargs):
    return 1 / catalog["radius"].to(u.arcmin)

def zweight(catalog, redshift_key, source_redshift):
    cat_redshifts = catalog[redshift_key]
    return cat_redshifts*source_redshift - cat_redshifts**2

