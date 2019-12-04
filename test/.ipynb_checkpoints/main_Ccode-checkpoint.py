import FilamentPaint
import sys
sys.path.insert(0,'../src')

import numpy as np
import matplotlib.pyplot as pl

from Functions import *
from Sky import Sky
from MagField import MagField
from FilPop import FilPop

output_tqumap	= 'test.fits'
nside			= 512
Npix_box		= 256
theta_LH_RMS	= None # in degrees
size_scale		= 0.7
size_ratio		= 0.25
Nfil			= 10000
size_box		= 1800.0 # physical size box
# Create the sky object
sky			= Sky(nside)
# Create the magnetic field object
magfield		= MagField(size_box,Npix_box,12345)
# Create the filament population object
population		= FilPop(Nfil,theta_LH_RMS,size_ratio,size_scale,magfield,fixed_distance=True)

for n in range(Nfil):
	TQUmap			= FilamentPaint.Paint_Filament(n,sky,population,magfield)
	tqumap			= np.array(TQUmap)
	sky.Tmap 	+= tqumap[0,:]
	sky.Qmap 	+= tqumap[1,:]
	sky.Umap 	+= tqumap[2,:]

sky.save_sky(output_tqumap)