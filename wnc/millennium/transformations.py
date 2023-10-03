from functools import partial
from loguru import logger
from typing import Callable
import astropy.units as u
from cosmap.analysis import CosmapBadSampleError
import json
from pathlib import Path
from dask.distributed import get_worker
from copy import copy
from itertools import product
import numpy as np

column_file = Path(__file__).parent / "columns.json"

with open(column_file, "r") as f:
    columns = json.load(f)


def build_filter(z_s, inner_radius, radius, limiting_magnitude, column_aliases) -> Callable:
    if "mag_i" not in column_aliases or "redshift" not in column_aliases:
        raise ValueError("Column aliases must contain 'mag_i' and 'redshift' keys")
    def filter_function(catalog, radius, limiting_magnitude):
        mask = catalog[column_aliases["mag_i"]] < limiting_magnitude
        mask & (catalog["radius"] > inner_radius)
        mask & (catalog["radius"] < radius)
        mask & (catalog[column_aliases["redshift"]] < z_s)
        return catalog[mask]
    return partial(filter_function, radius=radius, limiting_magnitude=limiting_magnitude)


class Setup:

    @staticmethod
    def build_filter_sets(redshifts, inner_radius, radius, limiting_magnitude) -> dict:
        dataset_name = "ms"
        column_aliases = columns[dataset_name]
        #make sure radius is a list
        if isinstance(radius, u.Quantity) and radius.isscalar:
            radius = [radius]
        #make sure limiting_magnitude is a list
        if isinstance(limiting_magnitude, float):
            limiting_magnitude = [limiting_magnitude]
        filters = {}
        #Create keys for the filter sets by combining the radius and limiting magnitude
        for name, redshift in redshifts.items():
            lens_filters = {}
            for r in radius:
                for m in limiting_magnitude:
                    #get limiting magnitude without any a period
                    mag_key = str(m).replace(".", "")
                    #drop any trailing zeros
                    mag_key = mag_key.rstrip("0")
                    #make a key by combing the radius and limiting magnitude
                    key = f"{int(r.value)}_{mag_key}"
                    filter = build_filter(redshift, inner_radius, r, m, column_aliases)
                    lens_filters[key] = filter
            filters[name] = lens_filters
        return filters
    
    @staticmethod
    def build_output_paths(filter_sets, base_output_path, fields_to_sample) -> dict:
        if fields_to_sample is not None:
            splits = [f.split("_") for f in fields_to_sample]
            field_tuples = [(int(s[0]), int(s[1])) for s in splits]
            print(fields_to_sample)
        else:
            print("NO FIELDS TO GET")
            exit()
            n = list(range(8))
            field_tuples = list(product(n, n))
        
        paths = {}
        for lens in filter_sets.keys():
            print(lens)
            lens_folder = base_output_path / lens
            if lens_folder.exists():
                logger.warning(f"Output path {lens_folder} already exists. Continuing anyway...")
            lens_folder.mkdir(parents=False, exist_ok=True)
            lens_filters = filter_sets[lens]
            folders = {f"{lens}_{f}": lens_folder / f for f in lens_filters.keys()}
            print(folders)
            for path in folders.values():
                if path.exists():
                    logger.warning(f"Output path {path} already exists. Continuing anyway...")
                path.mkdir(parents=False, exist_ok=True)
            for folder_key, folder_path in folders.items():
                keys = [f"{folder_key}_{i}_{j}" for i, j in field_tuples]
                file_paths = {k: folder_path / f"{k[-3:]}.csv" for k in keys}
                paths.update(file_paths)
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
    def apply_filters(catalog, filters, output_paths, sample_region, *args, **kwargs):
        distances = sample_region.coordinate.separation(catalog['coordinates'])
        catalog["radius"] = distances
        catalogs = {}
        for lens_name, lens_filters in filters.items():
            lens_catalogs = {}
            for filter_name, filter in lens_filters.items():
                filtered_catalog = filter(catalog)
                lens_catalogs.update({filter_name: filtered_catalog})
            catalogs.update({lens_name: lens_catalogs})
        return catalogs
    


    @staticmethod
    def count(catalogs, sample_region, current_field, redshifts, *args, **kwargs):
        output = {}
        current_field = get_worker().dataset.get_parameter("ms_field")
        for lens_name, lens_catalogs in catalogs.items():
            lens_z = redshifts[lens_name]
            for filter_name, catalog in lens_catalogs.items():
                catalog_outputs = compute_weightfns(catalog, "z_spec", lens_z)
                output_key = f"{lens_name}_{filter_name}_{current_field[0]}_{current_field[1]}"
                output[output_key] = {
                    "ra": sample_region.coordinate.ra.deg,
                    "dec": sample_region.coordinate.dec.deg,
                    **catalog_outputs
                }
        return output

                
        
    
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


def ncount(catalog, *args, **kwargs):
    return np.ones(len(catalog))

def oneoverr(catalog, *args, **kwargs):
    return 1 / catalog["radius"].to(u.arcmin)

def zweight(catalog, redshift_key, source_redshift):
    cat_redshifts = catalog[redshift_key]
    return cat_redshifts*source_redshift - cat_redshifts**2

