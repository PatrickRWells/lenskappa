from abc import ABC

from lenskappa.utils.datamanager import SurveyDataManager
from lenskappa.catalog import SkyCatalog2D
from lenskappa.spatial import SkyRegion
import logging

class DataSet(ABC):

    def __init__(self, name, *args, **kwargs):
        """
        An abstract base class for objects with multiple
        types of data. Example: A survey, which includes catalogs,
        bright star masks, and regions.  
        """
        self._name = name
        self._datamanager = SurveyDataManager(name)



class SkyDataSet(DataSet):

    def __init__(self, name, *args, **kwargs):
        super().__init__(name, *args, **kwargs)
        self._validated = False
        """
        Class for datasets on the sky (or datasets that can be projected onto the sky)
        Should not be directly instantiated.
        Can be counted
        """

    def _validate(self, *args, **kwargs):
        try:
            _cat = isinstance(self._catalog, SkyCatalog2D)
            if not _cat:
                logging.error("SkyDataSet requires a SkyCatalog2D")
                return
        except:
            logging.error("No catalog found!")
            return
    
        try:
            _reg = isinstance(self._region, SkyRegion)
            if not _reg:
                logging.error("SkyDataSet requires a SkyRegion")
                return
        except:
            logging.error("No region found!")
            return

        self._validated = True


