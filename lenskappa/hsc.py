import re
from astropy.coordinates import SkyCoord
import astropy.units as u
from shapely import geometry
import pandas as pd
import numpy as np
import math



class hsc_catalog:
    def __init__(self, file):
        self._data = pd.read_csv(file)
        print("Splitting data...")
        self.subsets = {id: self._data[self._data.tract == id] for id in self._data.tract.unique()}
        self._patch_data = {}

    
    def get_objects_within_patches(self, tracts):
        frames = []
        for id, data in tracts.items():
            patch_ints = [self._patch_tuple_to_int(key) for key in data['subregions'].keys()]
            for key in patch_ints:
                if key not in self._patch_data.keys():
                    mask = (self.subsets[id].patch).isin(patch_ints)
                    self._patch_data.update({key: self.subsets[id][mask]})
                frames.append(self._patch_data[key])
        return pd.concat(frames)

    def _patch_tuple_to_int(self, patch_tuple):
        if patch_tuple[0] == 0:
            return patch_tuple[1]
        else:
            return 100*patch_tuple[0] + patch_tuple[1]


class field:
    def __init__(self, top_region, sub_regions):
        self._region = top_region
        self._subregions = {key: region(value) for key,value in sub_regions.items()}
    def get_inside(self, center, aperture):
        #Given a circular aperture, determine which subregions it intersects with
        subregions = {}
        aperture_deg = aperture.to(u.deg).value
        circle_center = geometry.point.Point(center.ra.degree, center.dec.degree)
        aperture_poly = circle_center.buffer(aperture_deg)
        for id, region in self._subregions.items():
            inside = region.get_inside(aperture_poly)
            if inside or inside is None:
                subregions.update({id: inside})
        return subregions


    
    def tile(self, aperture_size = 120 * u.arcmin):
        #Tiles a region with the minimum number of aperture_size radius
        #circular sub-regions that covers the whole region
        #Note, this assumes a rectangular region
        bounds = self._region.bounds 
        size = aperture_size.to(u.degree).value
        dx = np.abs(bounds[0] - bounds[2])
        dy = np.abs(bounds[1] - bounds[3])
        num_x = math.ceil(dx/size) + 1
        num_y = math.ceil(dy/size) + 1
        bot_left_x, bot_left_y = min(bounds[0], bounds[2]), min(bounds[1], bounds[3])
        x_coords = [bot_left_x + i*size for i in range(num_x)]
        y_coords = [bot_left_y + i*size for i in range(num_y)]
        self._x_tiles = x_coords
        self._y_tiles = y_coords
        self._aperture_size = aperture_size
        
class region:
    def __init__(self, data):
        try: 
            self._type = data['type']
        except:
            self._type = 'None'

        self._center = SkyCoord(*data['center'], unit="deg")
        corner_ra, corner_dec = list(map(list, zip(*data['corners'])))
        self._corners = SkyCoord(corner_ra, corner_dec, unit="deg")
        if 'subregions' in data.keys():
            self._subregions = {key: region(value) for key, value in data['subregions'].items()}

    def get_inside(self, aperture_poly):
        subregions = {}
        if not hasattr(self, "_poly"):
            self._poly = geometry.Polygon(list(zip(self._corners.ra.degree, self._corners.dec.degree)))
        intersects = aperture_poly.intersects(self._poly)
        if not intersects:
            return False
        
        elif hasattr(self, "_subregions"):
            for id, subregion in self._subregions.items():
                inside = subregion.get_inside(aperture_poly)
                if inside or inside is None:
                    subregions.update({id: inside})
        return {'ref': self, 'subregions': subregions}
    
    def set_data(self, data):
        self._data = data


def read_hsc_tracts(tractfile):
    tracts = {}
    with open(tractfile) as tf:
        for line in tf:
            if line.startswith('*'):
                continue

            nums = re.findall(r"[-+]?\d*\.\d+|\d+", line)
            tract_num = int(nums[0])
            patch = re.search("Patch", line)
            if patch is None:
                if tract_num not in tracts.keys():
                    tracts.update({tract_num: {'corners': [], 'type': 'tract', 'subregions': {} } } )
                if 'Corner' in line:
                    tracts[tract_num]['corners'].append((float(nums[-2]), float(nums[-1])))
                elif 'Center' in line:
                    tracts[tract_num].update({'center': (float(nums[-2]), float(nums[-1]))})

            else:
                patch = re.findall(r'\d,\d', line)
                patch_val = tuple(map(int, patch[0].split(',')))
                if patch_val not in tracts[tract_num]['subregions'].keys():
                    tracts[tract_num]['subregions'].update({patch_val: {'corners': [], 'type': 'patch'}})
                if 'Corner' in line:
                    tracts[tract_num]['subregions'][patch_val]['corners'].append((float(nums[-2]), float(nums[-1])))
                elif 'Center' in line:
                    tracts[tract_num]['subregions'][patch_val].update({'center': (float(nums[-2]), float(nums[-1]))})
    return tracts
            
if __name__ == '__main__':    
    reg_points = [(40, -6.5), (29.5, -6.5), (29.5, -2), (40, -2)]
    sky_region = geometry.Polygon(reg_points)
    aperture = 120 * u.arcsec

    print("Loading tract data..")
    tracts = read_hsc_tracts('tracts_patches_W-w02.txt')
    sky = field(sky_region, tracts)
    print("Tiling sky")
    sky.tile(aperture)
    exit()

    center = SkyCoord(35.0573440373 , -4.55841069495, unit="deg")
    inside_tracts = sky.get_inside(center, aperture)
    print("Tiling sky...")

    print("Loading catalog data...")
    cat = hsc_catalog("/Users/patrick/Documents/Current/Research/LensEnv/HSC/catalogs/W2.csv")
    print("Mapping tracts and patches to field...")
