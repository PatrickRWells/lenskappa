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

def _compute(catfile, which='all', params={}):
    if which=='all':
        weights = toml.load('config/counts.toml')['counts']
    else:
        weights = which
    weight_fns = _load_weightfns(list(weights.keys())) #Load the weighting functions
    needspars = []
    for weight, pars in weights.items(): #Check to make sure the catalog contains all the required values
        if pars['params'] == 0:
            continue
        for par_val in pars['params']:
            if par_val.startswith('cat'):
                parname = par_val.split('.')[-1]
                print(parname)


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
    cfhtmask_file = "/Volumes/workspace/CFHT/masks/W1m0m0_izrgu_finalmask_mosaic.fits"
    cfhtcat_file = "/Volumes/workspace/CFHT/catalogs/W1m0m0/W1m0m0_24galphotmstar.cat"
    
    #compute_weights("HE0435", 'HE0435IRACbpz_nobeta_i24.cat', cfhtcat_file, cfhtmask_file, {'radius': 45})
    _compute('all')