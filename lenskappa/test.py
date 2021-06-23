
from astropy.coordinates import SkyCoord
import regions
from lenskappa.mask import RegMask


center = SkyCoord(45, 45, unit="deg")
new_center = SkyCoord(46, 46, unit="deg")
mask_file = "/Users/patrick/Documents/Current/Research/LensEnv/HSC/HSC-SSP_brightStarMask_Arcturus/reg/patches/10039/BrightStarMask-10039-0,0-HSC-I.reg"
mask = RegMask.from_file(mask_file,center)
new_mask = mask.relocate(new_center)
