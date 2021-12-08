from re import M
from lenskappa.catalog import SkyCatalog2D, params
from lenskappa.catalog import QuantCatalogParam, SingleValueParam
from lenskappa.datasets.surveys import hsc
from lenskappa.spatial import SkyRegion, CircularSkyRegion
from lenskappa.catalog import MaxValueFilter
from lenskappa.counting import RatioCounter


from astropy.coordinates import SkyCoord
import astropy.units as u
from shapely import geometry
import logging
import time

# Lenskappa expects standard names for catalog columns. It can handle
# Non-standard names using paramater objectes like these one.
# These can also handle other parameters associated with the analysis
# Like thr edshift of the source quasar

mass_param = QuantCatalogParam("demp_sm", "m_gal", is_log=True)
z_param = QuantCatalogParam("demp_photoz_best", "z_gal")
z_s_param = SingleValueParam("z_s", 1.523)

#Read in the catalog for your field of interest. Currently, only CSV files are supported
lens_field = SkyCatalog2D.read_csv("/home/prwells/data/J0924/lens_cat.csv", params=[mass_param, z_param, z_s_param])

# The HSC W02 only has full coverage in about half of its area
# So we create a box defining the edge of what we want to compare to
box = geometry.box(30, -6, 39, -2)
region = SkyRegion(box.centroid, box)
#Sky regions take the center of the region as an argument (I plan to remove this soon)

#Initialize the comparison field
survey = hsc("W02", params=[mass_param, z_param, z_s_param], frame=region)

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
magmax_filter = MaxValueFilter('i_cmodel_mag', 24)


#Initialize the counter object. mask=True means the counter will take the 
#bright star masks into account when computing the number counts
counter = RatioCounter(lens_field, survey, lens_region, mask=True)

#Add the catalog filters to the counter, which will automatically apply them.
#An "absolute" filter is applied to the catalogs before any other work is done.
#which = 'both' just means the filters will be to both the control and lens catalogs
counter.add_catalog_filter(zmax_filter, name='zmax')
counter.add_catalog_filter(magmax_filter, name='magmax')

counter.get_weight_ratios("all", output_file="lenskappa_test.csv", num_samples=1000, threads= 4, overwrite=True)
#The first parameter in get_weight_ratios indicates which weights to compute
#For now it's best to just leave this on all.
 
