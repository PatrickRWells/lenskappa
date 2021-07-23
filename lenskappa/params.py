
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

class QuantCatalogParam(CatalogParam):
    """
    Class for handling parameters with numerical data.
    Can deal with logs
    """
    def __init__(self, col_name, std_name, is_log = False, *args, **kwargs):
        super().__init__(self, col_name, std_name, *args, **kwargs)
        self._is_log = False

    def get_values(self, cat, *args, **kwargs):
        vals = super().get_values(cat, *args, **kwargs)
        if not self._is_log:
            return vals
        else:
            return np.power(10, vals)


    
