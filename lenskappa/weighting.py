# P Wells
# This code is based on earlier code by CE Rusu. It is used to calculate weights for objects in a lens field. 

from lenskappa.lens import Lens
import pandas as pd
import toml
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import time

class weight:
    def __init__(self, weightname, weight_config = None):
        self._name = weightname
        if weight_config is None:
            self._config = self._load_weight_config()
        else:
            self._config=weight_config
        self._load_weight()
        
    
    def _load_weight_config(self):
        import os
        loc = "config/counts.toml"
        path = os.path.join(os.path.dirname(os.path.abspath(lenskappa.__file__)), loc)
        return path
    
    def _load_weight(self):
        """
        Loads the appropriate weightfunction
        """
        from lenskappa import weightfns
        try:
            weightdata = self._config[self._name]
        except:
            print("Error: Weight {} not found in configuration files".format(self._name))
            return
        self._weightfn = getattr(weightfns, self._name)
        pars = weightdata['params']
        self._cat_params = [par.split('.')[-1] for par in pars if par.startswith('cat')]
        self._other_params = [par for par in pars if not par.startswith('cat')]
    
    def compute_weight(self, catalog, pars={}):
        """
        Computes the weights based on an input catalog
        parameters:
        catalog: Pandas dataframe or astropy table with catalog data
        pars: additional information, such as parameter maps
        """
        #Check to make sure the required parameters exist in the catalog
        #Or are mapped in pars
        cat_pars_found = [par in catalog.columns for par in self._cat_params]
        self._parmap = {}
        #Check for catalog params
        for parname, is_found in zip(self._cat_params, cat_pars_found):
            if is_found:
                self._parmap.update({parname: parname}) 
                continue
            try:
                parmap = pars['map'][parname]
                self._parmap.update({parname: parmap})
            except:
                print("Error: parameter {} required to calculate weight {} but couldn't find it".format(parname, self._name))
                return
        for parname in self._other_params:
            try:
                parval = pars[parname]
                self._parmap.update({parname: parval})
            except:
                print("Error: unable to find value for parameter {} required to calculate weight {}".format(parname, self._name))
                return
        masks = np.array([True]*len(catalog))
        for par_name, col_name in self._parmap.items():
            try:
                mask = ~pd.isna(catalog[col_name])
                masks = masks & mask
            except:
                pass

        return self._weightfn(catalog[masks], self._parmap)

def load_all_weights():
    import os
    import lenskappa
    #Loads all weights found in the given config file and returns a dictionary
    loc = "config/counts.toml"
    config = os.path.join(os.path.dirname(os.path.abspath(lenskappa.__file__)), loc)

    weight_config = toml.load(config)
    weights = {key: weight(key, weight_config) for key in weight_config.keys()}
    return(weights)

def get_weight_rations_hsc(lens_name, lens_catalog_file, hsc_catalog_file, hsc_file, aperture, field_params, hsc_params):
    pass



def get_weight_ratios_cfht(lens_name, lens_catalog_file, cfht_catalog_file, cfht_mask_file, aperture, field_params, cfht_params):
    from lenskappa.mask import CFHTMask
    from lenskappa.catalog import catalog
    from lenskappa.lens import Lens

    #Load the masks and catalogs
    cfht_mask = CFHTMask(cfht_mask_file)
    cfht_catalog = catalog(cfht_catalog_file)
    lens_catalog = catalog(lens_catalog_file, verbose=True)
    lens_obj = Lens(lens_name)

    #Get the distances between the objects in the lens catalgo and the lens
    #Setting objects closer than 10 arcseconds to 10
    lens_obj.get_distances(lens_catalog, min_limit = 10)

    outer_radius = aperture['outer']
    inner_radius = aperture['inner']

    lens_catalog.apply_aperture(inner_radius, outer_radius)

    try: weights = field_params['weights']
    except:
        print("Warning: weights not specified. Defaulting to all")
        weight = 'all'

    #Build the map of fields from the large scale CFHT mask
    print("Building fieldmap on CFHT mask")    
    cfht_mask.build_maskmap(outer_radius)

    #Get the pixel locations for the objects in CFHT catlog
    print("Building pixelmap for CFHT catalog")
    cfht_catalog.get_macro_pixelmap(cfht_mask)

    #Remove objects from the CFHT catalog that are masked
    print("Removing objects from CFHT catalog that are masked")
    #cfht_catalog.remove_masked(cfht_mask)


    if weight == 'all':
        weights = load_all_weights()
    views = cfht_mask.get_views()
    cfht_params.update({'internal': True, 'aperture': aperture})
    field_params.update({'internal': False, 'center': lens_obj['center_coords']})
    field_weights = {key: weight.compute_weight(lens_catalog, field_params) for key, weight in weights.items()}
    ratios = {key: [] for key in weights.keys()}
    num_views = len(views)
    skipped = 0
    start_time = time.time()
    for index in range(num_views):
        


        view = views[index]


        if index%100 == 0:
            print("Done {} out of {} views. Skipped {}".format(str(index), str(num_views), str(skipped)))
            print("TIME ELAPSED: {}".format(str(time.time() - start_time)))
        if not view.apply(max_coverage = 0.75):
            skipped += 1
            continue
        cfht_masked_cat = cfht_catalog.apply_mask(view, cfht_params)
        field_masked_cat = lens_catalog.apply_mask(view, field_params)
        cfht_weighted_counts = {key: weight.compute_weight(cfht_masked_cat, cfht_params) for key, weight in weights.items()}
        field_weighted_counts = {}
        for key, field_weight in field_weights.items():
            if type(field_weight) is pd.Series:
                field_weighted_counts.update({key: field_weights[key][field_masked_cat.index]})
            else:
                field_weighted_counts.update({key: weights[key].compute_weight(field_masked_cat, field_params)})


        temp_ratios = {}
        view.clear()
        missing = []
        for weight, cfht_count in cfht_weighted_counts.items():
            if type(cfht_count) is pd.Series:
                cfht_sum = sum(cfht_count)
                field_weighted_sum = sum(field_weighted_counts[weight])
                if cfht_sum == 0:
                    temp_ratios.update({weight: -1})
                else:
                    temp_ratios.update({weight: field_weighted_sum/cfht_sum})
            elif cfht_count is not None:
                if cfht_count == 0:
                    temp_ratios.update({weight: -1})
                else:
                    temp_ratios.update({weight: field_weighted_counts[weight]/cfht_count})
            elif cfht_count is None:
                weights.pop(weight)
                field_weights.pop(weight)
                ratios.pop(weight)
                missing.append(weight)
        for key, ratio in temp_ratios.items():
            if key in missing:
                continue
            ratios[key].append(ratio)



    return pd.DataFrame(ratios)
           