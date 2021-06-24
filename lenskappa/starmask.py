from abc import ABCMeta, abstractmethod
import logging
from astropy.coordinates import SkyCoord, SkyOffsetFrame
import regions
from shapely import geometry, affinity
import astropy.units as u
import pandas as pd
from copy import deepcopy

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
    def relocate(self, new_center):
        """
        Rotates the region catalog on the sky
        Returns a new Mask object with the set centered at the given value
        Parameters:
            to - SkyCoord marking the center of the ending location
        """
        pass

    @abstractmethod
    def from_file(cls, file, center):
        pass

    @abstractmethod
    def mask_catalog(self, catalog, region):
        pass
    
    @abstractmethod
    def mask_external_catalog(self, catalog, region):
        pass


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
        except:
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
        

    def relocate(self, new_center):
        """
        Recenters a particular bright star masks at a new location.
        Used primarily when comparing a control field to a lens field

        Parameter:
            new_center <SkyCoord>: New location to center the mask
        Returns:
            new_mask <RegStarMask>: Deepcopy of self, with center changed
        """

        sep = self._center.separation(new_center)
        pa = self._center.position_angle(new_center)
        new_mask = deepcopy(self)
        new_mask._shape_centers = new_mask._shape_centers.directional_offset_by(pa, sep)
        centers = new_mask._shape_centers
        for i, shape in enumerate(new_mask._shapely_regdata):
            shape.center.x = centers[i].ra.degree
            shape.center.y = centers[i].dec.degree
        new_mask._center = new_center
        
        return new_mask

        
class FitsStarMask(StarMask):
    pass
    # TODO


class StarMaskCollection(metaclass=ABCMeta):
    def __init__(self, masks):
        """
        Container class for collections of masks.
        Should contain methods that allow it to be used as if it was a mask
        """
        self._masks = masks
    
    @abstractmethod
    def mask_catalog(self, catalog, region):
        pass
    
    @abstractmethod
    def mask_external_catalog(self, catalog, region):
        pass
    
    