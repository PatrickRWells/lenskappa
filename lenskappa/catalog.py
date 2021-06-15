from astropy.coordinates.sky_coordinate import SkyCoord
import pandas as pd
import numpy as np
from astropy import wcs, coordinates, units
import astropy.units as u
from shapely import geometry
from copy import copy
import logging
from shapely.geometry import geo
import toml
from lenskappa import base_config
from abc import ABCMeta, abstractmethod




class Catalog(metaclass=ABCMeta):

    def __init__(self, cat, *args, **kwargs):
        self._cat = cat
        self._needed_parameter_types = {}
        self._parmap = {}
    
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
    
    def __str__(self):
        return self._cat.__str__()

    @staticmethod
    def require_validation(pars):
        """
        Decorator to ensure a particular set of parameters
        Have been validate against the input catalog
        eg:
        @Catalog.require_validation(['x', 'y'])
        def func(self, *args, **kwargs)

        This would ensure that the parameters x, and y have been mapped
        to the appropriate catalog columns and their types are valid

        """
        def wrapper(function):
            def check_validation(self, *args, **kwargs):
                try:
                    checkvalid = self._valid
                except:
                    logging.error("Catalog parameters have not been validated!")
                    return None
                
                isvalid = [checkvalid[name] for name in pars]
                if np.all(isvalid):
                    return function(self, *args, **kwargs)
                else:
                    logging.error("Catalog parameter(s) {} are not valid".format([par for index, par in enumerate(pars) if not isvalid[index]]))
        

            return check_validation
        return wrapper


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
            return cls(cat, *args, **kwargs)
        except:
            logging.error("Pandas could not read file {}".format(file))
            return None
    
    def load_params(self, input_map, *args, **kwargs):
        """
        Used for marking catalog columns with standard labels
        e.g. The code expects columns "ra" and "dec", while 
        the catalog is labeled with "RaMean" and "DecMean"

        Parameters:
            input_map {str: str} - key-value pairs where the key
                is the standard column name and the value is the
                actual column name.
        """



        try:
            parmap = self._parmap
        except:
            self._parmap = {}
        
        for par, parmap in self._parmap.items():
            try: 
                input_parmap = input_map[par]
                self._parmap.update({par: input_parmap})
            except:
                logging.warning("No rename found for parameter {}, defaulting to {}".format(par, parmap))
        self._validate_parmap()
    
    def _validate_parmap(self):
        """
        Checks to make sure that all necessary parameters are present
        In the catalog and the values are of the expected type.
        Note: Currently checks for type by attempting to cast the column as
        the expected type. This should handle problems like getting a
        np.float instead of float, but may have unintended side effects
        """
        self._valid = {}

        for par, partype in self._needed_parameter_types.items():

            single_parmap = self._parmap[par]
            try:
                col = self._cat[single_parmap]
            except:
                logging.warning("Unable to find column {} in data".format(single_parmap))
                self._valid.update({par: False})
                continue
            
            try:
                coldata = col.astype(partype)
                self._valid.update({par : True})
            except:
                logging.error("Unable to cast column {} as type {}".format(single_parmap, partype))
                self._valid.update({par : False})

            

    @abstractmethod
    def add_subregion(self, *args, **kwargs):
        pass

    @abstractmethod
    def filter_by_subregion(self, *args, **kwargs):
        pass


class Catalog2D(Catalog):
    def __init__(self, cat, *args, **kwargs):
        super().__init__(cat)
        self._needed_parameter_types.update({'x': float, 'y': float})
        self._points_initialized = False
        self._subregions = {}



    def add_subregion(self, name, polygon, *args, **kwargs):
        """
        Adds a subregion to the catalog
        This is mostly just a convinient way of subdividing catalog data

        params:
            name <str>: Name for the subregion
            polygon <shapely.geometry.Polygon>: Polygon defining the region

        """
        if type(polygon) == geometry.Polygon:
            self._subregions.update({name: {'poly': polygon, 'points': [] } } )
        else:
            logging.error("Tried to add a subregion, but the region was not a polygon")
            
        
    @Catalog.require_validation(['x', 'y'])
    def filter_by_subregion(self, name, *args, **kwargs):
        if name not in self._subregions.keys():
            logging.error("Subregion {} does not exist".format(name))
            return

        if not self._points_initialized:
            self._init_points()
        if not self._subregions[name]['points']:
            self._init_subregion(name)
    
    def _init_points(self):
        x = self._cat[self._parmap['x']]
        y = self._cat[self._parmap['y']]
        self._points = [geometry.Point(xy) for xy in list(zip(x,y))]

    def _init_subregion(self, name):
        pass

class SkyCatalog2D(Catalog2D):

    def __init__(self, cat, *args, **kwargs):
        super().__init__(cat)
        self._parmap.update({'x': 'ra', 'y': 'dec'})

    def _init_points(self):
        #Need to implement a version of this function that works with SkyCoords
        #Since Astropy has utilities that take care of 3D corrections 


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
    lens_coord = SkyCoord(141.23246, 2.32358, unit="deg")
    aperture = 120*u.arcsec

    reg_mask = geometry.Point(lens_coord.ra.degree, lens_coord.dec.degree)
    print(aperture.to(u.degree).value)
    reg_mask.buffer(aperture.to(u.deg).value)

    cat = SkyCatalog2D.read_csv('/Users/patrick/Documents/Current/Research/LensEnv/0924/weighting/lens_cat.csv')
    cat.load_params({'x': 'ra', 'y': 'dec'})
    cat.add_subregion("test", reg_mask)
    cat.filter_by_subregion("test")