# P Wells
# This code is based on earlier code by CE Rusu. It is used to calculate weights for objects in a lens field. 

import logging
import toml
import lenskappa

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
    
    def compute_weight(self, catalog, pars={}):
        """
        Computes the weights based on an input catalog
        parameters:
        catalog: Pandas dataframe or astropy table with catalog data
        pars: additional information, such as parameter maps
        """
        #Check to make sure the required parameters exist in the catalog
        #Or are mapped in pars
        cat_pars = catalog.get_parmap()
        #Check for catalog params
        for std_name in self._cat_params:
            if std_name not in cat_pars.keys():
                logging.error("Looked for parameter {} in the input catalog, but couldn't find it".format(std_name))
                exit()
        for parname in self._other_params:
            try:
                parval = cat_pars[parname]
            except:
                print("Error: unable to find value for parameter {} required to calculate weight {}".format(parname, self._name))
                exit()
        
        self._parmap = cat_pars

        return self._weightfn(catalog, self._parmap)    


def load_all_weights():
    import os
    import lenskappa
    #Loads all weights found in the given config file and returns a dictionary
    loc = "weighting/config/counts.toml"
    config = os.path.join(os.path.dirname(os.path.abspath(lenskappa.__file__)), loc)

    weight_config = toml.load(config)
    weights = {key: weight(key, weight_config) for key in weight_config.keys()}
    return(weights)

def load_some_weights(weight_names):
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
            weightfn = weight(name, weight_config)
            weights.update({name: weightfn})
        except:
            logging.warning("Unable to find weight function {}. Skipping...".format(name))
    if len(weights) != 0:
        return(weights)
    else:
        logging.error("No weights were loaded...")
        return None