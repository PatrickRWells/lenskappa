# P Wells
# This code is based on earlier code by CE Rusu. It is used to calculate weights for objects in a lens field. 

from lenskappa.lens import Lens
from lenskappa.mask import Mask
import toml
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u

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
        return self._weightfn(catalog, self._parmap)

def load_all_weights():
    import os
    import lenskappa
    #Loads all weights found in the given config file and returns a dictionary
    loc = "config/counts.toml"
    config = os.path.join(os.path.dirname(os.path.abspath(lenskappa.__file__)), loc)

    weight_config = toml.load(config)
    weights = {key: weight(key, weight_config) for key in weight_config.keys()}
    return(weights)                                                               
           
