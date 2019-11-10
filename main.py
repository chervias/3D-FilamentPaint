import sys
sys.path.insert(0,'src')
import numpy as np
import healpy as hp
import matplotlib.pyplot as pl
from Functions import *
from Sky import Sky
from MagField import MagField
from FilPop import FilPop
import multiprocessing as mp

output_tqumap	= 'test.fits'
nside			= 512
Npix_box		= 256
Nfil			= 7
theta_LH_RMS	= 10. # in degrees
size_box		= 1500.0 # physical size box
# Create the sky object
sky				= Sky(nside)
# Create the magnetic field object
magfield		= MagField(size_box,Npix_box,12345)
# Create the filament population object
population		= FilPop(Nfil,theta_LH_RMS,magfield)

data	= [(n,sky,population,magfield) for n in range(Nfil)]
p 		= mp.Pool(processes=2)
for x in p.imap_unordered(paint_filament_aux, data,chunksize=2):
	t,q,u		= x
	sky.Tmap 	+= t
	sky.Qmap 	+= q
	sky.Umap 	+= u

#sky.mask[pix_filament]	= 1.0
#hp.mollview(sky.mask,nest=False)
#hp.mollview(sky.Tmap,nest=False)
#hp.mollview(sky.Qmap,nest=False)
#hp.mollview(sky.Umap,nest=False)
#pl.show()

sky.save_sky(output_tqumap)
