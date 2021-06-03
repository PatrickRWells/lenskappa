from astropy.io import fits
from astropy.io.fits import column
import astropy.units as u
from astropy.wcs import wcs, utils
from astropy.coordinates import SkyCoord
import numpy as np
from copy import copy
import re
import pandas as pd
import regions

class mask:
    pass


class CFHTMask(mask):
    """
    Class for managing CFHT masks
    Includes logic to break mask into individual fields
    Each individual field is an object of type Mask that contains a view
    of the underlying cfht maskfile
    """
    def __init__(self, maskfile):
        self._maskfile = maskfile
        self._pix = 0.187*u.arcsec
    
    
    def __getitem__(self, index):
        return self._mask[0].data[index]

    
    def build_maskmap(self, radius):
        """
        Splits the underlying mask into a map of smaller masks, each with the given radius
        Adapted from code by CE Rusu
        """
        self._mask = fits.open(self._maskfile)
        self._wcs = wcs.WCS(self._mask[0].header)
        shape = self._mask[0].data.shape
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
        self._views = self.build_views(field_centers, radius, self._pix, self._wcs)
    
    def build_views(self, field_centers, radius, pix, wcs):
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
                        params = {'center': center, 'radius': rad_pix, 'pix': pix, 'mask': self._mask, 'wcs': wcs}
                        view = maskview(params)
                        views[i][j][k][l] = view
        return views
    
    def get_views(self):
        if not hasattr(self, '_views'):
            print("ERROR: Must build views before returning them")
        
        return self._views.flatten()
    
    def get_wcs(self):
        try:
            return self._wcs
        except:
            print("ERROR: Tried to get a WCS object from a mask but it has not been initialized!")
    def get_pixel_values(self, x_pix, y_pix):
        return self._mask[0].data[x_pix, y_pix]
   
    
class maskview:
    def __init__(self, params):
        self._view_type = 'circular' #Currently the only kind supported
        if self._view_type is 'circular':
            self._radius = params['radius']
            self._center = params['center']
            self._pix = params['pix']
            self._wcs = params['wcs']
        
        self._mask = params['mask']
        self._mask_header = self._mask[0].header


    def get_view_center(self):
        return utils.pixel_to_skycoord(*self._center, self._wcs), self._radius*self._pix
        
    def apply(self, max_coverage=0.75):

        if self._view_type == 'circular':
            return self._get_circular_view(max_coverage)
    
    def _get_circular_view(self, max_coverage=0.75):
        import time
        center_x = self._center[0]
        center_y = self._center[1]
        xlow, xhigh = int(round(center_x - self._radius)), int(round(center_x + self._radius))
        ylow, yhigh = int(round(center_y - self._radius)), int(round(center_y + self._radius))
        view = self._mask[0].data[xlow: xhigh, ylow: yhigh]
        flattened_view = view.flatten()
        self.maskfrac = len(flattened_view[flattened_view != 0])/len(flattened_view)
        if self.maskfrac > 0.75:
            return False

        idxs, idys = np.indices(view.shape)
        distxs, distys = idxs - view.shape[0]/2, idys - view.shape[1]/2
        dists = np.sqrt(distxs**2 + distys**2)
        inside = dists <= self._radius
        final_view = copy(view)

        final_view[~inside] = 8192
        import matplotlib.pyplot as plt
        #for idx, idy, val in ndennumerate
        self._view = final_view
        self.shape = final_view.shape

        return True
    def clear(self):
        del(self._view)
    def get_vals(self, xs, ys):
        return self._view[xs, ys]
    
    def get_macro_vals(self, xs, ys):
        return self._mask[0].data[xs, ys]

 
    def get_all(self):
        return self._view

    def is_flipped(self):
        return self._mask_header['CDELT1'] < 0


class HSCMask(mask):
    
    def __init__(self, file):
        if not file.endswith('.reg'):
            print("Error: expected a region file for a HSC Mask")
            return None
        labels = ['shape', 'ra', 'dec', 'ax1', 'ax2', 'angle']
        regions_list = []
        print("Reading in region file...")
        with open(file) as f:
            for line in f:
                if line.startswith('circle'):
                    series = self._parse_circle(line)
                    regions_list.append(series)
                if line.startswith('box'):
                    series = self._parse_box(line)
                    regions_list.append(series)
        self._regdata = pd.DataFrame(regions_list, columns = labels)
        self._coords = SkyCoord(self._regdata['ra'], self._regdata['dec'], unit="deg")
        print("done")
    def _parse_circle(self, line):
        shape = line.split("#")[0].strip()
        nums = [float(num) for num in re.findall(r'\d+\.*\d*', line)]
        center = (nums[0], nums[1])
        radius = nums[2]*u.degree
        return ["circle", *center, radius, radius, 0]



    def _parse_box(self, line):
        shape = line.split("#")[0].strip()
        nums = [float(num) for num in re.findall(r'\d+\.*\d*', line)]
        center = (nums[0], nums[1])
        r1, r2 = nums[2]*u.degree, nums[3]*u.deg
        angle = nums[4]
        if angle != 0:
            print("ANGLE: {}".format(angle))
        return ["box", *center, r1, r2, angle]

    def get_circular_view(self, center, radius=120*u.arcsec):
        distances = self._coords.separation(center).to(u.arcsec)
        in_view = self._coords.separation(center).to(u.arcsec) <= 2*radius
        pixscale = 0.5*u.arcsec
        num_pix = int(4*radius/pixscale)
        wcs_header = {
            'CTYPE1': 'RA---TAN',
            'CUNIT1': 'degree',
            'CDELT1': -pixscale.to(u.degree).value,
            'CRPIX1': num_pix/2,
            'CRVAL1': center.ra.degree,
            'NAXIS1': num_pix,
            'CTYPE2': 'DEC--TAN',
            'CUNIT2': 'degree',
            'CDELT2': pixscale.to(u.degree).value,
            'CRPIX2': num_pix/2,
            'CRVAL2': center.dec.degree,
            'NAXIS2': num_pix
        }
        wcs_temp = wcs.WCS(wcs_header)
        pix_array = np.ndarray(shape=(num_pix, num_pix))
        pix = self._build_pixelarray(in_view, wcs_temp, pix_array=pix_array, mask_center=center, mask_radius=radius)
        import matplotlib.pyplot as plt
        plt.imshow(pix, origin='lower', cmap='gray')
        plt.show()

    def _build_pixelarray(self,regmask, wcs_obj, pix_array, mask_center, mask_radius):
        for index, shape_line in self._regdata[regmask].iterrows():
            center_coord = SkyCoord(shape_line['ra'], shape_line['dec'], unit="deg")
            if shape_line['shape'] == 'circle':
                region_sky = regions.CircleSkyRegion(center_coord, shape_line['ax1'])
                mask = region_sky.to_pixel(wcs_obj).to_mask('center')
              
                pix_array += mask.to_image(pix_array.shape)
            if shape_line['shape'] == 'box':
                region_sky = regions.RectangleSkyRegion(center_coord, shape_line['ax1'], shape_line['ax2'], shape_line['angle']*u.degree)
                mask = region_sky.to_pixel(wcs_obj).to_mask('center')
                pix_array += mask.to_image(pix_array.shape)

        aperture = regions.CircleSkyRegion(mask_center, mask_radius)
        aperture_mask = aperture.to_pixel(wcs_obj).to_mask('center').to_image(pix_array.shape) == 0.
        pix_array[aperture_mask] = 1
        pix_array[pix_array != 0] = 8192
        shape = pix_array.shape
        left_bound, right_bound = int(pix_array.shape[0]/2 - pix_array.shape[0]/4), int(pix_array.shape[0]/2 + pix_array.shape[0]/4), 
        lower_bound, upper_bound = int(pix_array.shape[1]/2 - pix_array.shape[1]/4), int(pix_array.shape[1]/2 + pix_array.shape[1]/4) 

        pix_array = pix_array[left_bound:right_bound, lower_bound:upper_bound]
        return pix_array

if __name__ == "__main__":
    regfile = "/Users/patrick/Documents/Current/Research/LensEnv/HSC/HSC-SSP_brightStarMask_Arcturus/reg/tracts/BrightStarMask-9617-HSC-I.reg"
    center = SkyCoord(218.9928211, 0.7540274, unit="deg")
    hsc = HSCMask(regfile)
    hsc.get_circular_view(center, radius=180*u.arcsec)