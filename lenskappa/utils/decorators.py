import logging
import numpy as np

def require_points(fn):
    """
    Decorator for Catalog to check that the object points in a catalog have
    been initialized
    """
    def wrapper(self, *args, **kwargs):
        if not self._points_initialized:
            self._init_points()
        return fn(self, *args, **kwargs)
    return wrapper


def require_validation(params, *args, **kwargs):
    """
    Decorator for Catalog to ensure a particular set of parameters
    Have been validate against the input catalog
    eg:
    @require_validation(['x', 'y'])
    def func(self, *args, **kwargs)

    This would ensure that the parameters x, and y have been mapped
    to the appropriate catalog columns and their types are valid

    """
    def outer(fn, *args, **kwargs):

        def wrapper(self, *args, **kwargs):
            try:
                checkvalid = self._valid
            except:
                logging.warning("Catalog parameters have not been validated!")
                return None
            isvalid = [checkvalid[name] for name in params]
            if np.all(isvalid):
                return fn(self, *args, **kwargs)
            else:
                logging.error("Catalog parameter(s) {} are not valid".format([par for index, par in enumerate(params) if not isvalid[index]]))
        return wrapper
    return outer
