from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np
from copy import copy

class Mask:
    def __init__(self, maskfile, type = 'CFHT'):
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
            'CDELT2': self.pix.to(u.degree).value,
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

if __name__ == '__main__':
    from lens import Lens
    import pandas as pd
    cfhtmask_file = "/Volumes/workspace/CFHT/masks/W1m0m0_izrgu_finalmask_mosaic.fits"
    cfhtcat_file = "/Volumes/workspace/CFHT/catalogs/W1m0m0/W1m0m0_24galphotmstar.cat"
    lens = Lens("HE0435")
    mask = Mask(cfhtmask_file)
    cat = np.loadtxt('HE0435IRACbpz_nobeta_i24.cat')
    cat = pd.DataFrame(cat, columns=['uk1', 'uk2', 'ra', 'dec', 'mag1', 'magerr1', 'uk3', 'zphot', 'mag2', 'mag3'])
    masking = mask.apply_to_catalog(lens, cat, 90)
    print(masking[1][1][10][10])
