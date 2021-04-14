# P R Wells
# This code is based on earlier code by CE Rusu. It is used to calculate weights for objects in a lens field. 

from lens import Lens
import toml
import numpy as np
from astropy.table import Table
from astropy.io import fits


def load_cat(catfile):
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
def delete_masked(cat, mask):
    cat_mask = np.where(mask[0].data[cat['y_pix'].astype(int), cat['x_pix'].astype(int)] != 0)
    return cat[cat_mask]
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


if __name__ == "__main__":
    data = load_cat("HE0435IRACbpz_nobeta_i24.cat")
    mask = fits.open("mskHE0435_asecrad120.fits")
    delete_masked(data,mask)