# P Wells
# This code is based on earlier code by CE Rusu. It is used to calculate weights for objects in a lens field. 

from lens import Lens
from mask import Mask
import toml
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u

def load_lenscat(catfile):
    try:
        cat = np.loadtxt(catfile)
    except:
        raise IOError("File {} not found".format(catfile))
        exit()
    column_details = toml.load("config/catalog.toml")
    expected_num_cols = len(column_details['expected_cols'])
    found_num_cols = cat.shape[1]
    if expected_num_cols != found_num_cols:
        print("Warning: expected to find {} columns in catalog but found {}".format(expected_num_cols, found_num_cols))
        print("You may want to check your catalog against config/catalog.toml")
    else:
        print("Catalog loaded sucessfully.")
        print("Assumed catalog columns match with config/catalog.toml")
        return Table(data=cat, names=column_details['expected_cols'])

def load_cfhtcat(catfile):
    return catfile

def apply_cfht_mask(mask, cat, cells_on_a_side):
    pixCFHT = 0.187 * u.arcsec
    cfht_wcs = WCS(mask)
    if radius == 45: overlap = 2
    if radius == 60: overlap = 3
    if radius == 90: overlap = 4
    if radius == 120: overlap = 5
    centerfieldspix = np.zeros((overlap,overlap,2))
    unmaskedcell = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))

def get_halomass(cat):
    pass

def init_lens(name):
    return Lens(name)
def load_constants():
    config_location = "config/constants.toml"
    return toml.load(config_location)

def remove_from_cat():
    pass
    #Need to be able to process a catlog of objects to remove. 
def declare_counts(which='all'):
    pass



def compute_weights(lensid, lens_catfile, cfht_catfile, cfht_maskfile, params, lens_maskfile = None):
    lens = Lens(lensid) #Get the information for the lens
    lenscat = load_lenscat(lens_catfile) #Load the lens environment catalog
    cfhtcat = load_cfhtcat(cfht_catfile) #Load the CFHT catalog
    cfht_mask = Mask(cfht_maskfile) #Load the CFHT mask
    constants = load_constants() #Load constants
    lenscat = lens.get_distances(lenscat) #Get distances between lens and objects in its catalog
    radius = params['radius']
    catalog_masks = cfht_mask.apply_to_catalog(lens, lenscat, radius) #Get the catalog masks
    which = params.pop['which'] if 'which' in params.keys() else 'all'
    output = _compute_galweights(lens_catfile, catalog_masks, which, params)
    return output

def _compute_galweights(catfile, which='all', params={}):
    if which=='all':
        weights = toml.load('config/counts.toml')['counts']
    else:
        weights = which
    weight_fns = _load_weightfns(list(weights.keys())) #Load the weighting functions
    allpars = []
    weight_needspars = {}
    for weight, pars in weights.items(): #Get the catalog parameters that are needed to caclulate the weights
        if pars['params'] == 0:
            continue
        weightpars = []
        for par_val in pars['params']:
            if par_val.startswith('cat'):
                parname = par_val.split('.')[-1]
                if parname not in allpars:
                    allpars.append(parname)
                weightpars.append(parname)
        weight_needspars.update({weight: weightpars})        
    #Check to see if the catalog parameters are either in the catalog, or mapped in the params
    par_assignment = {}
    for par_name in allpars:
        if par_name in catfile.columns:
            par_assignment.update({par_name: par_name})
        else:
            try:
                if par_name in params['par_assign']:
                    par_assignment.update({par_name: params['par_assign'][par_name]})
                else: print("Error: unable to find parameter {} in catalog".format(par_name))
            except:
                print("Error: unable to find parameter {} in catalog".format(par_name))
    weights_final = []
    for weight, parvals in weight_needspars.items():
        if all([par in par_assignment.keys() for par in parvals]):
            weights_final.append(weight)
    output = {}
    for weights_key in weights_final:
        print("Calculating weights for {}".format(weights_key))
        weight_vals = weight_fns[weights_key](catfile, par_assignment, params['other_data'])
        output.update({weights_key: weight_vals})
    return output

                
        


def _load_weightfns(weights):
    import weightfns
    fns = {}
    for key in weights:
        try:
            fns.update({key: getattr(weightfns, key)})
        except:
            print("Error: Unable to find weight function \'{}\' in weightfns.py. Skipping...".format(key))
    return fns
                                                               
           


if __name__ == "__main__":
    import pandas as pd
    from copy import copy
    cfhtmask_file = "/Volumes/workspace/CFHT/masks/W1m0m0_izrgu_finalmask_mosaic.fits"
    cfhtcat_file = "/Volumes/workspace/CFHT/catalogs/W1m0m0/W1m0m0_24galphotmstar.cat"
    cat = pd.read_csv('/Volumes/workspace/LensENV/SDSS0924/photometry/hsc/187280.csv')
    assign = {'z_gal': 'photoz_best', 'm_gal': 'stellar_mass', 'r': 'dist'}
    lens = Lens("J0924")
    print("GETTING DISTANCES")
    cat = lens.get_distances(cat)
    print("GOT DISTANCES")
    cat = copy(cat[(cat['photoz_best'] <= lens['z_s']) & (cat['i_cmodel_mag'] <= 24.0)])
    cat.to_csv('catalog.csv')
    other_data = {'z_s': lens['z_s']}
    pars = {'par_assign': assign, 'other_data': other_data}
    weights = _compute_galweights(cat, 'all', pars)
    weight_frame = pd.DataFrame()
    for weight_key, weights in weights.items():
        print(type(weights))
        if type(weights) == np.float64:
            print("PASS")
            pass
        else:
            weight_frame[weight_key] = weights
    weight_frame.to_csv("weights.csv")
