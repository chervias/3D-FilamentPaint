import FilamentPaint
import sys
sys.path.insert(0,'code')

import numpy as np
import matplotlib.pyplot as pl

from Sky import *
from MagField import MagField
from FilPop import FilPop

output_tqumap	= 'test.fits'
nside			= 2048
Npix_box		= 256
theta_LH_RMS	= -1 # in degrees
size_scale		= 0.1
size_ratio		= 0.25
Nfil			= 100
size_box		= 400.0 # physical size box
slope			= 2.4

r_unit_vectors  = get_r_unit_vectors(nside)
local_triad     = get_local_triad(nside,r_unit_vectors)

# Create the magnetic field object
magfield		= MagField(size_box,Npix_box,12345)
# Create the filament population object
population		= FilPop(Nfil,theta_LH_RMS,size_ratio,size_scale,slope,magfield,1234,fixed_distance=True)

TQUmap_total = np.zeros((3,12*nside**2))
for n in range(Nfil):
	TQUmap			= FilamentPaint.Paint_Filament(n,nside,local_triad,population,magfield)
	for i in range(3):
		TQUmap_total[i,:] 	+= TQUmap[i,:]

hp.write_map(output_tqumap,TQUmap_total,nest=False,overwrite=True)