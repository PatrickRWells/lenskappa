from lenskappa import Counter, SkyCatalog2D
from lenskappa.surveys import hsc
from lenskappa.region import SkyRegion, CircularSkyRegion
from lenskappa.filter import MaxValueFilter


from astropy.coordinates import SkyCoord
import astropy.units as u
from shapely import geometry
import logging
import time

# Lenskappa expects standard names for catalog columns. It can handle
# Non-standard names using paramater maps like this one.
# The key is the standard name, while the value is the name in the catalog
# Paramater maps also contain parameters that are needed to compute weights,
# Such as the redshift of the source (z_s)
lens_parmap = {'m_gal': 'demp_sm', 'z_gal': 'demp_photoz_best', 'z_s': 1.523}
field_parmap = {'m_gal': 'demp_sm', 'z_gal': 'demp_photoz_best', 'z_s': 1.523}

m_gal_

#Read in the catalog for your field of interest. Currently, only CSV files are supported
lens_field = SkyCatalog2D.read_csv("/home/prwells/data/J0924/lens_cat.csv", parmap=lens_parmap)

# The HSC W02 only has full coverage in about half of its area
# So we create a box defining the edge of what we want to compare to
box = geometry.box(30, -6, 39, -2)
region = SkyRegion(box.centroid, box)
#Sky regions take the center of the region as an argument (I plan to remove this soon)

#Initialize the comparison field
survey = hsc("W02", parmap=field_parmap, frame=region)

#Define the region corresponding to the lens field
#The radius of this region will be used when comparing to the control field
aperture = 120*u.arcsec
center = SkyCoord(141.23246, 2.32358, unit="deg")
lens_region = CircularSkyRegion(center, aperture)

#Place a filter on the catalogs
#The first parameter is the name in the actual catalog, or the standard name (assuming it has been mapped as above)
#This example removes all objects with i-band magnitude > 24
#And all objects more distance than z = 1.523
zmax_filter = MaxValueFilter('z_gal', 1.523)
magmax_filter = MaxValueFilter('i_cmodel_mag', 1.523)


#Initialize the counter object. mask=True means the counter will take the
#bright star masks into account when computing the number counts
counter = Counter(lens_field, survey, lens_region, mask=True)

#Add the catalog filters to the counter, which will automatically apply them.
#An "absolute" filter is applied to the catalogs before any other work is done.
#which = 'both' just means the filters will be to both the control and lens catalogs
counter.add_catalog_filter(zmax_filter, name='zmax', filter_type='absolute', which='both')
counter.add_catalog_filter(magmax_filter, name='magmax', filter_type='absolute', which='both')

counter.get_weight_ratios("all", output_file="lenskappa_test.csv", num_samples=1000, threads= 4, overwrite=True)
#The first parameter in get_weight_ratios indicates which weights to compute
#For now it's best to just leave this on all.
 
