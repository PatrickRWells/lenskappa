
from abc import ABC, abstractmethod
import logging
import numpy as np

class Param(ABC):
    def __init__(self, *args, **kwargs):
        pass

class CatalogParam(Param):

    def __init__(self, col_name: str, std_name: str,  *args, **kwargs):
        """
        A class for handling catalog parameters. Note, this particular class DOES NOT
        check that the data frame actually contains the column until the values are requested.

        Arguments:
        col_name <str>: Name of the column in the dataframe
        std_name <str>: Standard name of the column

        """
        super().__init__(*args, **kwargs)
        self._col_name = col_name
        self._std_name = std_name
    
    def get_values(self, cat, *args, **kwargs):
        try:
            return cat[self._col_name]

        except:
            logging.error("Unable to find values for paramater {} in catalog".format(self._std_name))
            raise

    @property
    def standard(self):
        return self._std_name
    @property
    def col(self):
        return self._col_name

class QuantCatalogParam(CatalogParam):
    """
    Class for handling parameters with numerical data.
    Can deal with logs
    Can also accept an astropy unit.
    The reason for this is that Pandas dataframes has some weird buggy
    interactions with astropy units.
    """
    def __init__(self, col_name, std_name, unit = None, is_log = False, *args, **kwargs):
        super().__init__(col_name, std_name, *args, **kwargs)
        self._is_log = False
        self._unit = unit
    
    @property
    def unit(self):
        return self._unit
    def get_values(self, cat, *args, **kwargs):
        vals = np.array(super().get_values(cat, *args, **kwargs))
        if self._is_log:
            vals = np.power(10, vals)
        if self._unit is not None:
            return vals*self._unit
        else:
            return vals

class SingleValueParam(CatalogParam):
    
    def __init__(self, name, value):
        super().__init__(name, name)
        self._value = value

    def get_values(self, *args, **kwargs):
        return self._value
