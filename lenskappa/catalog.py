import pandas as pd
import numpy as np
from astropy import wcs, coordinates, units
import astropy.units as u
from copy import copy

class catalog:
    
    def __init__(self, cat, delimeter = ','):
        if type(cat) == str:
            self._data = pd.read_csv(cat, delimeter)
        if type(cat) == pd.DataFrame:
            self._data = cat
        self.columns = self._data.columns
        self._params = []
        self._pixmap_initialized = False
        self._masked_removed = False
        
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
            is_masked = pix_values != 0
            num_removed = len(self._data[is_masked])

            self._data.drop(self._data[is_masked].index, inplace=True)
            self._coords  = coordinates.SkyCoord(self._data['ra'], self._data['dec'], unit=(units.deg, units.deg))
            self._masked_removed = True
            print("Dropped {} out of {} objects from catalog, based on the mask".format(str(num_removed), str(num)))
        
        else:
            print("ERRROR: Tried to apply a mask to a catalog but pixel positions of catalog objects have not been initialized!")
    
    def apply_mask(self, mask, params, plot=False):
        if self._masked_removed and self._pixmap_initialized and params['internal']:
            return self._get_region(mask, params, plot)
        

        from copy import copy
        if not 'pix_x' in self.columns:
            self._init_pixelmap(mask, params)

        shape = self._wcs.array_shape

        x_rounded, y_rounded = round(self._data['pix_x']).astype(int), round(self._data['pix_y']).astype(int)

        in_frame = (x_rounded >= 0) & (y_rounded >= 0) & (x_rounded < shape[0] - 1) & (y_rounded < shape[1] - 1)
        return_cat = copy(self._data[in_frame])
        not_masked = mask.get_vals(x_rounded[in_frame], y_rounded[in_frame]) == 0
        if plot:
            self._plot_mask_withcat(mask, return_cat, not_masked)
        return return_cat[not_masked]
    
    def _get_region(self, mask, params, plot):
        center, radius = mask.get_view_center()
        data_copy = copy(self._data)
        data_copy['dist'] = self._coords.separation(center).to(u.arcsec)
        data_copy.drop(data_copy[data_copy['dist'] >= radius].index, inplace=True)
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

        