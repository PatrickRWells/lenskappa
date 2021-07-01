import logging
import regions
import astropy.units as u
import pandas as pd
import numpy as np

from shapely import geometry, affinity
from astropy.coordinates import SkyCoord
from abc import ABCMeta, abstractmethod

class StarMask(metaclass=ABCMeta):

    def __init__(self, data, center):
        """
        Class for defining masks
        This class should be used for a mask covering a relatively small region of sky
        For larger sets, use a MaskCollection

        """
        self._data = data
        self._center = center
    

    @abstractmethod
    def from_file(cls, file, center):
        pass

    @abstractmethod
    def mask_catalog(self, catalog, region):
        """
        Use when masking a catalog that falls INSIDE this bright star mask
        """
        pass

    def get_bool_region_mask(self, catalog, region):
        """
        Same as above, except it returns the boolean mask instead of the catalog itself.
        
        """
    
    @abstractmethod
    def mask_external_catalog(self, catalog, region):
        """
        Use when masking a catalog that falls outside this bright star mask.
        Mask will be rotated onto the catalog
        """


class RegStarMask(StarMask):

    def __init__(self, data, center):
        """
        Class for masks which are defined using region files.
        Originally built for use with HSC Survey masks
        """
        super().__init__(data, center)
        
        self._convert_data_to_shapely()
    

    @classmethod
    def from_file(cls, file, center):
        import regions
        try:
            regdata = regions.read_ds9(file)
        except Exception as e:
            print(e)
            logging.error("Unable to read maskfile {}\n"\
                          "Expected a region file".format(file))
        obj = cls(regdata, center)
        obj._file = file
        return obj

    def _convert_data_to_shapely(self):
        """
        Converts data originally read in from a .reg file to a shapely polygon.
        Shapely is much more feature-complete and well optimized than the astropy regions package
        And while technically it only support operations on a cartesian plane, this is not a problem
        for bright star masks that are a few arcseconds in size 
        Currently supports circles and rectangles
        """
        new_regions = []
        centers = []

        for regdata in self._data:

            regtype = type(regdata)
            centers.append((regdata.center.ra.degree, regdata.center.dec.degree))

            if regtype == regions.CircleSkyRegion:
                new_regions.append(self._parse_shapely_circle(regdata))
    
            elif regtype == regions.RectangleSkyRegion:
                new_regions.append(self._parse_shapely_rectangle(regdata))

            else:
                logging.error("Unexpected region type {} found. This region will be skipped".format(regtype))
                centers.pop(-1)
        self._shapely_regdata = pd.Series(new_regions)
        self._shape_centers = SkyCoord(centers, unit="deg")
        #TODO: Some sort of validation

    def _parse_shapely_circle(self, regdata):
        # Shapely circles are not true circles, but they are close
        # enough that it shouldn't cause any issues
        x_  = regdata.center.ra.degree
        y_ = regdata.center.dec.degree
        radius_ = regdata.radius.to(u.deg).value
        shape = geometry.Point(x_, y_)
        shape = shape.buffer(radius_)
        return shape

    def _parse_shapely_rectangle(self, regdata):
        x = regdata.center.ra.degree
        y = regdata.center.dec.degree
        width_ = regdata.width.to(u.deg).value
        height_ = regdata.height.to(u.deg).value
        angle_ = regdata.angle.to(u.deg).value
        x_min, x_max = x - width_/2, x + width_/2
        y_min, y_max = y - height_/2, y + height_/2
        shape = geometry.box(x_min, y_min, x_max, y_max)
        if angle_ != 0:
            shape = affinity.rotate(shape, angle_)
        return shape

    def get_bool_region_mask(self, catalog, region):
        """
        Given a region in the mask, and a catalog in the same region, return a 
        boolean mask for the catalog, where True indicates the object falls behind a mask
        and False indicates the object does not.
        
        """
        cat_points = catalog.get_points()
        masks_in_region = [reg for reg in self._shapely_regdata if region.intersects(reg)]
        mask = np.array([False]*len(catalog))
        for index, point in enumerate(cat_points):
            #Loop through the objects in the catalog
            if not mask[index]:
                #If they have not already been found to fall behind a mask
                for submask in masks_in_region:
                    #Check to see if they fall behind any of the masks
                    if submask.contains(point):
                        mask[index] = True
                        break
                        #If we find they are masked, move on to the next one
        return mask

    def mask_catalog(self, catalog, region):
        pass

    def mask_external_catalog(self, catalog, region):
        pass

        
class FitsStarMask(StarMask):
    pass
    # TODO


class StarMaskCollection:

    def __init__(self, masks):
        """
        Container class for collections of masks.
        Should contain methods that allow it to be used as if it was a mask
        """
        self._masks = masks


    def mask_catalog(self, catalog, region, subregions = None):
        input_cat = catalog
        if subregions is not None:
            for subregion in subregions:
                input_cat = self._masks[subregion].mask_catalog(input_cat)
        
        pass
    
    def mask_external_catalog(self, catalog, region):
        pass

    