from abc import ABCMeta, abstractmethod
import logging

class Region(metaclass=ABCMeta):
    def __init__(self, center, region, *args, **kwargs):
        self._center = center
        self._region = region
    
    def add_subregion(self, name, center, input_region, override = False, *args, **kwargs):
        """
        Adds a subreegion to the current region.
        This class also supports adding sub-subregions through the "propogate_subregion" method
        Subregions must be fuly enclosed by their parent regions. This can be overriden, but do so at your
        own peril
        Parameters:
            name: Name for the subregion. Can be anything that works as a dictionary key
            center: Tuple with (x,y) coordinates of the center
            input_region: Polygon or list of (x,y) points defining the corners
            override: Whether to override the interior requirement (see ablove)
        Returns:
            True if the subregion was sucessfully added, flase otherwise
        """
        try:
            subregions = self._subregions
        except:
            self._subregion_centers = {}
            self._subregions = {}

        within = self._region.contains(input_region)
        if not within and not override:
            logging.error("Tried to initialize subregion {}, but it is not contained with its parent region".format(name))
            return False
        elif isinstance(input_region, Region):
            self._subregions.update({name: input_region})
        else:
            self._subregions.update({name: self.build_region(center, input_region, *args, **kwargs)})
        self._subregion_centers.update({name: center})
        return True
    
            
    @staticmethod
    def build_region(center, input_region, *args, **kwargs):
        pass

    @abstractmethod
    def overlaps(self, second_region, *args, **kwargs):
        pass

class SkyRegion(Region):
    def __init__(self, center, region, *args, **kwargs):
        super().__init__(center, region, *args, **kwargs)

    def build_region(self, center, input_region, *args, **kwargs):
        return SkyRegion(center, input)

    def overlaps(self, second_region):
        return False

