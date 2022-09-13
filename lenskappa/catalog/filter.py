from abc import ABC, abstractmethod
import logging
import numpy as np


class Filter(ABC):

    def __init__(self, filter_type=None, *args, **kwargs):
        """
        Filters are used to modify the contents of a catalog.
        An example usage would be to remove objects beyond a certain redshift
        Or remove objects during weighting that fall too far from the center
        of the region.
        """
        self._type = filter_type
    
    @abstractmethod
    def __call__(self, catalog, *args, **kwargs):
        """
        Filters should implement a __call__ method.
        In order to apply the filter, call
            filter(catalog)
        Params:
        catalog <catalog.Catalog>: The catalog to be filtered

        Returns:
        filtered_catalog <catalog.Catalog>: The filtered catalog        
        """
        pass

class ColumnFilter(Filter):

    def __init__(self, column, *args, **kwargs):
        """
        A filter operating on a single column in a catalog.
        """
        super().__init__("column", *args, **kwargs)
        self._column = column
    
    @abstractmethod
    def __call__(self, catalog, *args, **kwargs):
        pass
    
class ColumnLimitFilter(ColumnFilter):

    def __init__(self, column, min=None, max=None, *args, **kwargs):
        """
        Place a limit on the value of a column numerical column
        Rows where the given column is outside the limit will be removed

        Params:
            column: The column to be filtered by
            min: the minimum allowed value
            max: the maximum allowed value
        """

        if min is None and max is None:
            logging.error("A column limit filter must have either a minimum or maximum value!")
            return
        super().__init__(column, *args, **kwargs)
        self._min = min
        self._max = max
    
    def __call__(self, catalog, *args, **kwargs):
        column = catalog[self._column]
        filter = np.ones(len(column), dtype=bool)
        if self._min is not None:
            filter = filter & (column > self._min)
        if self._max is not None:
            filter = filter & (column < self._max)
        
        if np.all(filter):
            return catalog
        else:
            return catalog[filter]
    

class MaxValueFilter(ColumnLimitFilter):
    
    def __init__(self, column, max, *args, **kwargs):
        """
        Convinience class for when only a maximum value is necessary
        """

        super().__init__(column, max = max, *args, **kwargs)
    
class MinValueFilter(ColumnLimitFilter):
    
    def __init__(self, column, min, *args, **kwargs):
        """
        Convinience class for when only a minimum value is necessary
        """

        super().__init__(column, min = min, *args, **kwargs)

class ColumnLimitFilterWithReplacement(ColumnLimitFilter):

    def __init__(self, column, min = None, max = None, *args, **kwargs):
        """
        This class is identical to a column limit filter
        Except that instead of removing objects that fall
        outside the bounds, it replaces their value with the
        value of the bound.
        Everything larger than max will be replaced with max
        And everything smaller than min will be replaced with min 
        """
        super().__init__(column, min, max)
    
    def __call__(self, catalog):
        new_catalog = catalog
        col = catalog[self._column]

        if self._min is not None:
            min_filter = (col > self._min)
            new_catalog = new_catalog.replace_values_by_mask(min_filter, self._column, self._min)
        
        if self._max is not None:
            max_filter = (col < self._max)
            new_catalog = new_catalog.replace_values_by_mask(max_filter, self._column, self._max)
        
        return new_catalog