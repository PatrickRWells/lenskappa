# P Wells
# This code is based on earlier code by CE Rusu. It is used to calculate weights for objects in a lens field. 

import logging
import toml
import lenskappa
import numpy as np

class weight:

    def __init__(self, weightname, weight_config = None, params = []):
        self._name = weightname
        self._params = params
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
        from lenskappa.weighting import weightfns
        try:
            weightdata = self._config[self._name]
        except:
            print("Error: Weight {} not found in configuration files".format(self._name))
            return
        self._weightfn = getattr(weightfns, self._name)
        pars = weightdata['params']
        self._cat_params = [par.split('.')[-1] for par in pars if par.startswith('cat')]
        self._other_params = [par for par in pars if not par.startswith('cat')]

        try:
            self._post = self._handle_post_fn(weightdata['post'])
        except:
            self._post = None
    
    def _handle_post_fn(self, name):
        if not name.startswith("np"):
            logging.error("Unable to load postprocess function for weight {}".format(self._name))
        else:
            funcname = name.split(".")[-1]
            try:
                return getattr(np, funcname)
            except:
                logging.error("Unable to load postprocess function for weight {}".format(self._name))
                return None

    def compute_weight(self, catalog, meds=False):
        """
        Computes the weights based on an input catalog
        parameters:
        catalog: Pandas dataframe or astropy table with catalog data
        pars: additional information, such as parameter maps
        """
        #Check to make sure the required parameters exist in the catalog
        #Or are mapped in pars
        #Check for catalog params
        for std_name in self._cat_params:
            catalog[std_name]
        for parname in self._other_params:
            try:
                parval = self._params[parname]
            except:
                print("Error: unable to find value for parameter {} required to calculate weight {}".format(parname, self._name))
                exit()
        weights = self._weightfn(catalog, **self._params)
        if meds and len(weights) != 0: #Second condition avoids crash
            weights = np.median(weights)*np.ones(len(weights), dtype=np.float64)
        if self._post is not None:
            weights = self._post(np.sum(weights))
        return weights
    
    def _compute_weight_sampled_params(self, catalog):
        """
        Compute weights when one or more
        params has a sample instead of a single value
        I realized too late that this will not work :P
        """
        samples = {key: catalog.get_samples(key) for key in self._sampled}
        if len(samples.keys()) != 1:
            logging.error("Currently unable to handle more than one sampled parameter")
            return self._weightfn(catalog)
        
        for name, sample_obj in samples.values():
            num_samples = sample_obj.num_samples
            storage = np.zeros()
            actual_vals = catalog[name]
            for sample in sample_obj.get_samples():
                catalog[name] = sample


            catalog[name] = actual_vals

        


def load_all_weights(params = []):
    import os
    import lenskappa
    #Loads all weights found in the given config file and returns a dictionary
    loc = "weighting/config/counts.toml"
    config = os.path.join(os.path.dirname(os.path.abspath(lenskappa.__file__)), loc)

    weight_config = toml.load(config)
    weights = {key: weight(key, weight_config, params) for key in weight_config.keys()}
    return(weights)

def load_some_weights(weight_names, params = []):
    import os
    import lenskappa
    #Loads all weights found in the given config file and returns a dictionary
    loc = "config/counts.toml"
    config = os.path.join(os.path.dirname(os.path.abspath(lenskappa.__file__)), loc)
    weight_config = toml.load(config)
    not_found = []
    weights = {}
    for name in weight_names:
        try:
            weightfn = weight(name, weight_config, params)
            weights.update({name: weightfn})
        except Exception as e:
            print(e) 
            logging.warning("Unable to find weight function {}. Skipping...".format(name))
    if len(weights) != 0:
        return(weights)
    else:
        logging.error("No weights were loaded...")
        return None

if __name__ == '__main__':
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from lenskappa.catalog import QuantCatalogParam, SingleValueParam, SkyCatalog2D

    ws = load_all_weights()
    m_gal_param = QuantCatalogParam('demp_sm', 'm_gal', is_log=True)
    z_gal_param = QuantCatalogParam('demp_photoz_best', 'z_gal')
    z_s_param = SingleValueParam('z_s', 1.523)
    ra_param = QuantCatalogParam('ra', 'ra', u.deg)
    dec_param = QuantCatalogParam('dec', 'dec', u.deg)

    pars = [m_gal_param, z_gal_param, z_s_param, ra_param, dec_param]

    lens_field = SkyCatalog2D.read_csv("/Users/patrick/Documents/Current/Research/LensEnv/0924/weighting/lens_cat.csv", params=pars)
    center = SkyCoord(141.23246, 2.32358, unit="deg")
    others = SkyCoord(lens_field['ra'], lens_field['dec'], unit="deg")
    distances = center.separation(others).to(u.arcsec).value
    r_param = QuantCatalogParam("dist", "r", unit=u.arcsec)
    lens_field['dist'] = distances
    lens_field.add_param(r_param)

    
    for name, w in ws.items():
        w1 = w.compute_weight(lens_field)
        w1_meds = w.compute_weight(lens_field, meds=True)
        print("Weight {}: {},   {}".format(name, w1, w1_meds))