import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from copy import copy

class catalog:
    
    def __init__(self, cat, delimeter = ','):
        if type(cat) == str:
            self._data = pd.read_csv(cat, delimeter)
        if type(cat) == pd.DataFrame:
            self._data = cat
        self.columns = self._data.columns
        
    def __str__(self):
        return self._data.__str__()
    def __getitem__(self, key):
        return self._data[key]
    def __len__(self):
        return(len(self._data))

    def __setitem__(self, key, value):
        self._data[key] = value
    
    def apply_mask(self, mask, params, plot=False):
        if params['internal']:
            return self._apply_mask_internal(mask, params, plot)
        

        from copy import copy
        if not 'x_pix' in self.columns:
            self._init_pixelmap(mask, params)

        shape = self._wcs.array_shape
        x_rounded, y_rounded = round(self._data['x_pix']).astype(int), round(self._data['y_pix']).astype(int)

        in_frame = (x_rounded >= 0) & (y_rounded >= 0) & (x_rounded < shape[0]) & (y_rounded < shape[1])
        return_cat = copy(self._data[in_frame])
        not_masked = mask.get_vals(x_rounded[in_frame], y_rounded[in_frame]) == 0
        if plot:
            self._plot_mask_withcat(mask, return_cat, not_masked)
        return return_cat[not_masked]

    def _apply_mask_internal(self, mask, params, plot=False):
        #Applies a mask to a catalog of the same region
        center, radius = mask.get_view_center()
        coords = SkyCoord(self._data['ra'], self._data['dec'], unit=(u.deg, u.deg))
        data_copy = copy(self._data)
        data_copy['dist'] = coords.separation(center).to(u.arcsec)
        data_copy.drop(data_copy[data_copy['dist'] >= radius].index, inplace=True)
        params['center'] = center
        x_pix, y_pix = self._return_pixelmap(data_copy, mask, params)
        x_rounded, y_rounded = round(x_pix).astype(int), round(y_pix).astype(int)
        not_masked = mask.get_vals(x_rounded, y_rounded) == 0
        if plot:
            self._plot_mask_withcat(mask, data_copy, not_masked)
        return data_copy[not_masked]

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
        self._data['x_pix'] = pix_x
        self._data['y_pix'] = pix_y

    def _return_pixelmap(self, cat, mask, params):
        tempcat = catalog(cat)
        tempcat._init_pixelmap(mask, params)
        return tempcat['x_pix'], tempcat['y_pix']

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

        plt.scatter(cat['x_pix'], cat['y_pix'], c=c, s=5)
        plt.show()

        