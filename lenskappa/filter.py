from abc import ABCMeta, abstractmethod
import logging
import numpy as np


class Filter:

    def __init__(self, filter_type=None, *args, **kwargs):
        """
        Filters are used to modify the contents of a catalog.
        An example usage would be to remove objects beyond a certain redshift
        Or remove objects during weighting that fall too far from the center
        of the file.
        """
        self._type = filter_type
    
    @abstractmethod
    def __call__(self, catalog, parmap, *args, **kwargs):
        """
        Filters shoudl implement a __call__ method.
        In order to apply the filter, call
            filter(catalog, parmap)
        Params:
        catalog <catalog.Catalog>: The catalog to be filtered
        parmap <dict>: The parameter map for the catalog.
            Filters by default look for standard column names.        
        
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
    def __call__(self, catalog, parmap, *args, **kwargs):
        pass
    
    def _get_column_values(self, catalog, parmap):
        try:
            colname = self._check_column_name(catalog, parmap)
            return catalog[colname]
        except:
            raise

        
    def _check_column_name(self, catalog, parmap):
        try:
            col = catalog[self._column]
            return self._column
        except:
            try:
                colname = parmap[self._column]
                col = catalog[colname]
                return colname
            except:
                logging.error("Unable to locate column {}".format(self._column))
                raise

class ColumnLimitFilter(ColumnFilter):

    def __init__(self, column, min=None, max=None, *args, **kwargs):
        """
        Place a limit on the value of a column.
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
    
    def __call__(self, catalog, parmap, *args, **kwargs):

        column = self._get_column_values(catalog, parmap)
        filter = np.ones(len(column), dtype=bool)
        if self._min is not None:
            filter = filter & (column > self._min)
        if self._max is not None:
            filter = filter & (column < self._max)
        
        if np.all(filter):
            return catalog
        else:
            return catalog.apply_boolean_mask(filter)
    

class MaxValueFilter(ColumnLimitFilter):
    
    def __init__(self, column, max, *args, **kwargs):

        super().__init__(column, max = max, *args, **kwargs)
    
class MinValueFilter(ColumnLimitFilter):
    
    def __init__(self, column, min, *args, **kwargs):

        super().__init__(column, min = min, *args, **kwargs)

class ColumnLimitFilterWithReplacement(ColumnLimitFilter):

    def __init__(self, column, min = None, max = None, *args, **kwargs):
        """
        This class is identical to a column limit filter
        Except that instead of removing objects that fall
        outside the bounds, it replaces their value with the
        value of the bound.
        Everything larger than min will be replaced with max
        And everything smaller than min will be replaced with min 
        """
        super().__init__(column, min, max)
    
    def __call__(self, catalog, parmap):
        new_catalog = catalog
        col = self._get_column_values(new_catalog, parmap)
        colname = self._check_column_name(new_catalog, parmap)

        if self._min is not None:
            min_filter = (col > self._min)
            new_catalog = new_catalog.replace_values_by_mask(min_filter, self._column, self._min)
        
        if self._max is not None:
            max_filter = (col < self._max)
            new_catalog = new_catalog.replace_values_by_mask(max_filter, self._column, self._max)
        
        return new_catalog


if __name__ == "__main__":
    from lenskappa import SkyCatalog2D

    cat = SkyCatalog2D.read_csv("/Users/patrick/Documents/Current/Research/LensEnv/0924/weighting/lens_cat.csv")
    filter = ColumnLimitFilterWithReplacement('i_cmodel_mag', min=20, max=24)
    parmap = {}

    new_cat = filter(cat, parmap)
    import matplotlib.pyplot as plt
    plt.hist(cat['i_cmodel_mag'], 50)
    plt.show()
    plt.hist(new_cat['i_cmodel_mag'],50)
    plt.show()
    print(len(cat))
    print(len(new_cat))



