import toml
import re
import astropy.units as u
import astropy.coordinates as coords



class Lens(dict):
    def __init__(self, name):
        self._config = "config/lens_data.toml"
        self._load_data(name)
    def _load_data(self, name):
        config_data = toml.load(self._config)
        if name not in config_data.keys():
            raise ValueError("Data for lens {} not found in configuration files".format(name))
        else:
            self.update({'name': name})
            self.update(config_data[name])
        self._parse_coords()
        self._parse_units()
    def _parse_coords(self):
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
        unit_items = [val for val in self.keys() if re.search(r'unit', val)]
        for unit in unit_items:
            key_val = unit.replace('_unit', '')
            unit_str = self[unit]
            unit_obj = getattr(u, unit_str)
            self[key_val] *= unit_obj
            self.pop(unit)

if __name__ == "__main__":
    data = Lens("HE1104")
    print(data)
