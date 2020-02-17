import FilamentPaint
import sys
sys.path.insert(0,'code')

import numpy as np
import matplotlib.pyplot as pl

from Functions import *
from Sky import Sky
from MagField import MagField
from FilPop import FilPop

output_tqumap	= 'test.fits'
nside			= 512
Npix_box		= 256
theta_LH_RMS	= -1 # in degrees
size_scale		= 10.0
size_ratio		= 0.25
Nfil			= 100
size_box		= 400.0 # physical size box
slope			= 2.4
# Create the sky object
sky			= Sky(nside)
# Create the magnetic field object
magfield		= MagField(size_box,Npix_box,12345)
# Create the filament population object
population		= FilPop(nside,Nfil,theta_LH_RMS,size_ratio,size_scale,slope,magfield,1234,fixed_distance=True)

for n in range(Nfil):
	print(n,population.sizes[n,2])

exit()

for n in range(Nfil):
	TQUmap			= FilamentPaint.Paint_Filament(n,sky,population,magfield)
	sky.Tmap 	+= TQUmap[0,:]
	sky.Qmap 	+= TQUmap[1,:]
	sky.Umap 	+= TQUmap[2,:]

sky.save_sky(output_tqumap)
