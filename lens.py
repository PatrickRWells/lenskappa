import toml
import re
import astropy.units as u
import astropy.coordinates as coords



class Lens(dict):
    def __init__(self, name, lensdata=None):
        """
        Basic class representing the lens system.
        Can be constructed with Lens("lens_name"), where lens name is one of the lenses
        in your database.
        """
        if lensdata is None:
            self._config = self._find_config() #"config/lens_data.toml"
        else:
            self._config = lensdata()
        self._load_data(name)
    def _find_config(self):
        """
        Internal function. Finds the location of the lens database.
        """
        import os
        import lenskappa
        loc = "config/lens_data.toml"
        path = os.path.join(os.path.dirname(os.path.abspath(lenskappa.__file__)), loc)
        return path



    def _load_data(self, name):
        """Internal function, used for loading lens data"""
        config_data = toml.load(self._config)
        if name not in config_data.keys():
            raise ValueError("Data for lens {} not found in configuration files".format(name))
        else:
            self.update({'name': name})
            self.update(config_data[name])
        self._parse_coords()
        self._parse_units()
        print("Data for lens {} loaded sucessfully".format(name))

    def _parse_coords(self):
        """Parses coordinates in lens configuration files and returns them as
           SkyCoord objects"""
        #Find coordinates and construct SkyCoord objects
        coord_keys = [val for val in self.keys() if re.search(r'coords',val) and not re.search(r'unit', val)]
        unit_keys = ['_'.join([coord, 'units']) for coord in coord_keys]
        for index, coord in enumerate(coord_keys):
            if not unit_keys[index] in self.keys():
                raise ValueError("No units provided for coordinate {}".format(coord))
            else:
                units = self[unit_keys[index]]
                unit_obj = [getattr(u, unit_str) for unit_str in units]
                coord_val = self[coord]
                frame_key = '_'.join([coord, 'frame'])
                if frame_key in self.keys():
                    frame = self[frame_key]
                    self.pop(frame_key)
                else:
                    frame = 'fk5'

                coord_obj = coords.SkyCoord(coord_val, unit=unit_obj, frame=frame)
                self[coord] = coord_obj
                self.pop(unit_keys[index])

    def _parse_units(self):
        """
        Internal function.
        Parases non-coordinate units found in lens configuration files and returns the appropriate
        object from astropy.units
        """
        unit_items = [val for val in self.keys() if re.search(r'unit', val)]
        for unit in unit_items:
            key_val = unit.replace('_unit', '')
            unit_str = self[unit]
            unit_obj = getattr(u, unit_str)
            self[key_val] *= unit_obj
            self.pop(unit)

    def get_distances(self, cat, ra_unit=u.degree, dec_unit=u.degree, append_to_cat = True):
        """
        Gets the distance in arcseconds between the lens and all objects in a given catalog.
        These can be appended to the catalog or returned seperately.
        """
        field_coords = coords.SkyCoord(cat['ra'], cat['dec'], unit=(ra_unit, dec_unit), frame='fk5')
        field_distances = field_coords.separation(self['center_coords'])
        if append_to_cat:
            cat['dist'] = field_distances.to(u.arcsec)
            return cat
        else:
            return field_distances.arcsec

    def get_location(self):
        return
