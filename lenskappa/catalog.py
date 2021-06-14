import pandas as pd
import numpy as np
from astropy import wcs, coordinates, units
import astropy.units as u
from shapely import geometry
from copy import copy
import logging
import toml
from lenskappa import base_config

class Catalog:

    def __init__(self, cat, param_map = {}, *args, **kwargs):
        self._cat = cat
        self._param_map = param_map
        #self._config = base_config.get_submodule_config(Catalog)
        self._init_cat(*args, **kwargs)
    
    def __getitem__(self, key):
        """
        Allows for masking as with a usual dataframe
        """
        return self._cat[key]

    def __setitem__(self, key, data):
        """
        Allows for setting values as with a usual dataframe
        """
        self._cat[key] = data
    
    def __len__(self):
        """
        Returns the length of the underlying data
        """
        return len(self._cat)
    
    @classmethod
    def read_csv(cls, file, config = None, *args, **kwargs):
        """Read in a file from a CSV
        Arguments:
            - file (str): path to CSV file
            - config (dict): dictionary of config parameters
        Rerturns:
            catalog.Catalog
        """
        logging.info("Reading in file {}".format(file))
        try:
            cat = pd.read_csv(file)
            return cls(cat, config, *args, **kwargs)
        except:
            logging.error("Pandas could not read file {}".format(file))
            return None
    

    def _init_cat(self, *args, **kwargs):
        pass

    def add_subregion(self, subregion_name, subregion_polygon, *args, **kwargs):
        if not hasattr(self, "_subregions"):
            self._subregions = {}
        self._subregions.update({subregion_name: subregion_polygon})  

    def filter_catalog_by_subregion(self, subregion_name):
        try: catpoints = self._catpoints
        except: self._init_catpoints()
    
    def _init_catpoints(self):
        ras = self._cat[self._param_map['ra']].to(u.deg)
        decs = self._cat[self._param_map['dec']].to(u.deg)
        


class catalog_old:
    
    def __init__(self, cat):
        if cat is not pd.DataFrame:
            print("Error: catalog expects a dataframe object")
            return
        self._data = cat
    
    def frame(self, frame):
        if not hasattr(self, "points"):
            self.init_points()
    
    def init_points(self):
        x = self._data.ra
        y = self._data.dec
        self._points = [geometry.Point(x_i, y[idx]) for idx, x_i in x]


        
    def __str__(self):
        return self._data.__str__()
    def __getitem__(self, key):
        return self._data[key]
    def __len__(self):
        return(len(self._data))

    def __setitem__(self, key, value):
        self._data[key] = value

    def remove_masked(self, mask):
        if self._pixmap_initialized:
            num = len(self._data)

            pix_values = mask.get_pixel_values(self._data['pix_x'], self._data['pix_y'])
            import matplotlib.pyplot as plt
            plt.hist(pix_values, [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60])
            plt.show()
            is_masked = (pix_values != 0)
            from collections import Counter
            print(Counter(is_masked))
            num_removed = len(self._data[is_masked])

            self._data.drop(self._data[is_masked].index, inplace=True)
            self._coords  = coordinates.SkyCoord(self._data['ra'], self._data['dec'], unit=(units.deg, units.deg))
            self._masked_removed = True
            print("Dropped {} out of {} objects from catalog, based on the mask".format(str(num_removed), str(num)))
        
        else:
            print("ERRROR: Tried to apply a mask to a catalog but pixel positions of catalog objects have not been initialized!")
    
    def apply_mask(self, mask, params, plot=False):
        if self._pixmap_initialized and params['internal']:
            return self._get_region(mask, params, plot)
    
        from copy import copy

        if (not hasattr(self, '_wcs')) or (not 'pix_x' in self.columns):
            self._init_pixelmap(mask, params)

        shape = self._wcs.array_shape

        x_rounded, y_rounded = round(self._data['pix_x']).astype(int), round(self._data['pix_y']).astype(int)
        in_frame = (x_rounded >= 0) & (y_rounded >= 0) & (x_rounded < shape[0] - 1) & (y_rounded < shape[1] - 1)
        return_cat = copy(self._data[in_frame])
        not_masked = mask.get_vals(x_rounded[in_frame], y_rounded[in_frame]) == 0
        if plot:
            print("PLOTTING 2")
            self._plot_mask_withcat(mask, return_cat, not_masked)
        return return_cat[not_masked]
    
    def apply_aperture(self, inner_radius, outer_radius):
        if 'dist' not in self._data.keys():
            print("ERROR: Trying to remove objects based on distances, but catalog does not have them")
            return
        
        mask = (self._data['dist'] < outer_radius) & (self._data['dist'] > inner_radius)
        self._data.drop(self._data[~mask].index, inplace=True)
    
    def _get_region(self, mask, params, plot):
        center, radius = mask.get_view_center()
        inner_radius = params['aperture']['inner']
        data_copy = copy(self._data)
        data_copy['dist'] = self._coords.separation(center).to(u.arcsec)
        mask = (data_copy['dist'] <= radius) & (data_copy['dist'] >= inner_radius) 
        data_copy.drop(data_copy[~mask].index, inplace=True)
        if len(data_copy) == 0:
            print(center)
        return data_copy

    def get_macro_pixelmap(self, mask):
        self._wcs = mask.get_wcs()
        self._coords  = coordinates.SkyCoord(self._data['ra'], self._data['dec'], unit=(units.deg, units.deg))
        pix_x, pix_y = self._wcs.world_to_pixel(self._coords)
        self._data['pix_x'] = np.round(pix_x).astype(int)
        self._data['pix_y'] = np.round(pix_y).astype(int)
        self._params.append('macro')
        self._pixmap_initialized = True

    def _init_pixelmap(self, mask, params):
        from astropy import wcs, coordinates, units
        if 'center' not in params.keys():
            print("Error: Need to generate pixel values but don't know where to put the center")
            return

        cdelt1, cdelt2 = mask._mask_header['CDELT1'], mask._mask_header['CDELT2']
        self._center = params['center']
        width, height = mask.shape
        wcs_input = {
            'CTYPE1': 'RA---TAN',
            'CUNIT1': 'deg',
            'CDELT1': np.abs(cdelt1),
            'CRPIX1': width/2,
            'CRVAL1': self._center.ra.degree,
            'NAXIS1': width,
            'CTYPE2': 'DEC--TAN',
            'CUNIT2': 'deg',
            'CDELT2': np.abs(cdelt2),
            'CRPIX2': height/2,
            'CRVAL2': self._center.dec.degree,
            'NAXIS2': height
        }
        self._wcs = wcs.WCS(wcs_input)
        self._coords = coordinates.SkyCoord(self._data['ra'], self._data['dec'], unit=(units.deg, units.deg))
        pix_x, pix_y = self._wcs.world_to_pixel(self._coords)
        self._data['pix_x'] = pix_x
        self._data['pix_y'] = pix_y

    def _return_pixelmap(self, cat, mask, params):
        tempcat = catalog(cat)
        tempcat._init_pixelmap(mask, params)
        return tempcat['pix_x'], tempcat['pix_y']

    def _plot_mask_withcat(self, mask, cat, catmask):
        import matplotlib.pyplot as plt
        maskdata = mask.get_all()
        if mask.is_flipped():
            maskdata = maskdata.T

        plt.imshow(maskdata, cmap='gray', origin='lower')
        c = ['blue']*len(cat)
        for index, masked in enumerate(catmask):
            if not masked:
                c[index] = 'red'

        plt.scatter(cat['pix_x'], cat['pix_y'], c=c, s=5)
        plt.show()



if __name__ == "__main__":
    pd.read_csv
    cat = Catalog.read_csv("/Users/patrick/Documents/Current/Research/LensEnv/0924/weighting/lens_cat.csv")
    mask = [False]*len(cat)
