from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np
from copy import copy

class Mask:
    def __init__(self, maskfile, view = None):
        self._mask = fits.open(maskfile)
        self.type = type
        if self.type == 'CFHT':
            self.pix = 0.187*u.arcsec
    def __getitem__(self, index):
        return self._mask[index]
    
    def construct_equivalent_pix_catalog(self, lens, catalog, radius):
        '''
        Construct equivalent pixel positions for a catalog that did not come from a particular image.
        Used when masking out objects based on the control field masks
        The lens is placed at the center of the so-called "image"
        '''
        lens_coords = lens['center_coords']
        width = height = int((2*radius/self.pix).value)
        wcs_input = {
            'CTYPE1': 'RA---TAN',
            'CUNIT1': 'deg',
            'CDELT1': self.pix.to(u.degree).value,
            'CRPIX1': width/2,
            'CRVAL1': lens_coords.ra.degree,
            'NAXIS1': width,
            'CTYPE2': 'DEC--TAN',
            'CUNIT2': 'deg',
            'CDELT2': -1.0*self.pix.to(u.degree).value,
            'CRPIX2': height/2,
            'CRVAL2': lens_coords.dec.degree,
            'NAXIS2': height
        }
        wcs_new = WCS(wcs_input)
        coords = SkyCoord(catalog['ra'], catalog['dec'], unit=(u.deg, u.deg))
        pix_x, pix_y = wcs_new.world_to_pixel(coords)
        catalog['pix_x'], catalog['pix_y'] = pix_x, pix_y
        return(catalog)
    
    def apply_to_catalog(self, lens, catalog, radius, zeropoint=0):
        if self.type == 'CFHT':
            return self._apply_to_catalog_cfht(lens, catalog, radius, zeropoint)
    def _apply_to_catalog_cfht(self, lens, catalog, radius, zeropoint=0):
        """
        Applies a given mask to a catalog, removing objects that fall within a masked area.
        Note: This code was written to break up very large masks from CFHTLens into smaller regions
        A future version may do the breaking up elsewhere.

        """

        pix = self.pix
        cells_on_a_side = int((len(self._mask[0].data[1]) * pix.value) / (2 * radius)) - 1
        if radius == 45: overlap = 2
        if radius == 60: overlap = 3
        if radius == 90: overlap = 4
        if radius == 120: overlap = 5
        if not 'dist' in catalog.columns:
            catalog = lens.get_distances(catalog)
        in_range = catalog['dist'] <= radius
        centerfieldspix = np.zeros((overlap,overlap,2))
        unmaskedcell = np.zeros((overlap,overlap,cells_on_a_side,cells_on_a_side))

        if 'pix_x' not in catalog.columns:
            catalog = self.construct_equivalent_pix_catalog(lens, catalog, radius)
        masking = np.ndarray(shape=(overlap, overlap, cells_on_a_side, cells_on_a_side, len(catalog)), dtype=bool )
        for k in range(overlap):
            for l in range(overlap):
                for i in range(cells_on_a_side):
                    for j in range(cells_on_a_side):
                        #Define the center of each cell as a matrix of pixel centers and choose which cells to discard because of large masked area

                        xlow = (2.0 * radius * (u.arcsec / pix).value) * i + (2.0 * radius * (u.arcsec / pix).value) * k / overlap   # coordinates of cell edges in the field mask; here x axis refers to the array axis, so the image y axis
                        xhigh = (2.0 * radius * (u.arcsec / pix).value) * i + (2.0 * radius * (u.arcsec / pix).value) + (2.0 * radius * (u.arcsec / pix).value) * k / overlap
                        ylow = (2.0 * radius * (u.arcsec / pix).value) * j + (2.0 * radius * (u.arcsec / pix).value) * l / overlap
                        yhigh = (2.0 * radius * (u.arcsec / pix).value) * j + (2.0 * radius * (u.arcsec / pix).value) + (2.0 * radius * (u.arcsec / pix).value) * l / overlap
                        centerfieldspix[k][l][0] = xlow + (xhigh - xlow) / 2.0 # coordinates of cell centers in the field mask, in CFTLENS-size pixels
                        centerfieldspix[k][l][1] = ylow + (yhigh - ylow) / 2.0

                        x_lenscat, y_lenscat = catalog['pix_x'], catalog['pix_y'] # copy the pixel coordinates
                        lenscoords_fieldx = float(xlow)+(1.0 * y_lenscat/((2/pix.value) * radius))*(float(xhigh)-float(xlow)) # project the lens catalogue onto the field mask; this is good, because it matches the formula for unmaskedfieldx
                        lenscoords_fieldy = float(ylow)+(1.0 * x_lenscat/((2/pix.value) * radius))*(float(yhigh)-float(ylow))
                        lenscoords_field = self._mask[0].data[lenscoords_fieldx.astype(int),lenscoords_fieldy.astype(int)]
                        if np.all(lenscoords_field != 0):
                            continue
                        mask = (lenscoords_field == 0) & (in_range) #Create a mask of the catalog that removes
                        masking[k][l][i][j] = mask                  #Objects that are either out of range or covered by a field mask
                        #lens_field_mask = np. # remove objects inside a field mask
                        #lensbpz_masked = lensbpz_masked[:-1] # delete the last column

        return masking
    
class maskview:
    def __init__(self, params):
        self._view_type = 'circular' #Currently the only kind supported
        if self._view_type is 'circular':
            self._radius = params['radius']
            self._center = params['center']
            self._pix = params['pix']
    
    def get_view_center(self):
        from astropy.wcs import WCS, utils
        temp_wcs = WCS(self._mask_header)
        return utils.pixel_to_skycoord(*self._center, temp_wcs), self._radius*self._pix
        
    def apply_to_mask(self, mask):
        self._mask_header = mask[0].header
        if self._view_type == 'circular':
            return self._get_circular_view(mask)
    
    def _get_circular_view(self, mask):
        center_x = self._center[0]
        center_y = self._center[1]
        xlow, xhigh = int(round(center_x - self._radius)), int(round(center_x + self._radius))
        ylow, yhigh = int(round(center_y - self._radius)), int(round(center_y + self._radius))
        view = mask[0].data[xlow: xhigh, ylow: yhigh]
        idxs, idys = np.indices(view.shape)
        distxs, distys = idxs - view.shape[0]/2, idys - view.shape[1]/2
        dists = np.sqrt(distxs**2 + distys**2)
        inside = dists <= self._radius
        final_view = np.zeros(view.shape)
        for (i,j), val in np.ndenumerate(view):
            final_view[i,j] = val if inside[i,j] else 8192
        #for idx, idy, val in ndennumerate
        self._view = final_view
        self.shape = final_view.shape
    def get_vals(self, xs, ys):
        return self._view[xs, ys]
    def get_all(self):
        return self._view

    def is_flipped(self):
        return self._mask_header['CDELT1'] < 0

class CFHTMask:
    """
    Class for managing CFHT masks
    Includes logic to break mask into individual fields
    Each individual field is an object of type Mask that contains a view
    of the underlying cfht maskfile
    """
    def __init__(self, maskfile):
        self._maskfile = maskfile
        self._pix = 0.187*u.arcsec
    
    def build_maskmap(self, radius):
        """
        Splits the underlying mask into a map of smaller masks, each with the given radius
        Adapted from code by CE Rusu
        """
        self._mask = fits.open(self._maskfile)
        shape = self._mask[0].data.shape
        print(shape)
        cells_on_a_side = int((len(self._mask[0].data[1]) *self._pix.value) / (2 * radius)) - 1
        if radius not in [45, 60, 90, 120]:
            print("Error, radius should be one of [45, 60, 90, 120]")
            return
        if radius == 45: overlap = 2
        if radius == 60: overlap = 3
        if radius == 90: overlap = 4
        if radius == 120: overlap = 5
        field_centers = np.zeros(shape=(cells_on_a_side, cells_on_a_side, overlap, overlap, 2))
        for k in range(overlap):
            for l in range(overlap):
                for i in range(cells_on_a_side):
                    for j in range(cells_on_a_side):
                        #Define the center of each cell as a matrix of pixel centers and choose which cells to discard because of large masked area
                        pix_diameter = 2.0*radius*(u.arcsec/self._pix).value
                        xlow = (pix_diameter) * i + (pix_diameter) * k / overlap   # coordinates of cell edges in the field mask; here x axis refers to the array axis, so the image y axis
                        xhigh = (pix_diameter) * (i + 1) + (pix_diameter) * k / overlap
                        ylow = (pix_diameter) * j + (pix_diameter) * l / overlap
                        yhigh = (pix_diameter) * ( j + 1)  + (pix_diameter) * l / overlap
                        field_centers[i][j][k][l][0] = xlow + (xhigh - xlow) / 2.0 # coordinates of cell centers in the field mask, in CFTLENS-size pixels
                        field_centers[i][j][k][l][1] = ylow + (yhigh - ylow) / 2.0
        print("Cells built")
        self._views = self.build_views(field_centers, radius, self._pix)
    
    def build_views(self, field_centers, radius, pix):
        shape = list(field_centers.shape)
        shape.pop(-1)
        rad_pix = radius/pix.value
        if len(shape) != 4:
            print("Error: Expected 4D input for fields")
        views = np.ndarray(shape=shape, dtype=maskview)

        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    for l in range (shape[3]):
                        center = field_centers[i][j][k][l]
                        params = {'center': center, 'radius': rad_pix, 'pix': pix}
                        view = maskview(params)
                        views[i][j][k][l] = view
        return views





if __name__ == '__main__':
    from lens import Lens
    import pandas as pd
    import matplotlib.pyplot as plt
    cfhtmask_file = "/Volumes/workspace/CFHT/masks/W1m0m0_izrgu_finalmask_mosaic.fits"
    cfhtcat_file = "/Volumes/workspace/CFHT/catalogs/W1m0m0/W1m0m0_24galphotmstar.cat"
    lens = Lens("HE0435")
    mask = CFHTMask(cfhtmask_file)
    mask.build_maskmap(radius=45)
    view = mask._views[10][10][1][1].apply_to_mask(mask._mask)
    #image = plt.imshow(view, cmap='gray')
    plt.show()
    #masking = mask.apply_to_catalog(lens, cat, 90)
