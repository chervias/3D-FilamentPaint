import FilamentPaint
import sys
sys.path.insert(0,'code')

import numpy as np
import matplotlib.pyplot as pl

from Sky import *
from MagField import MagField
#from FilPopCopy import FilPop
from FilPop import FilPop

output_tqumap	= 'test.fits'
dust_template = '/home/chervias/CMB-work/Filaments/COM_CompMap_IQU-thermaldust-gnilc-unires_2048_R3.00.fits'
nside			= 2048
Npix_box		= 256
theta_LH_RMS	= 55 # in degrees
size_scale		= 0.05
size_ratio		= 0.25
Nfil			= 2
size_box		= 400.0 # physical size box
slope			= 2.4

r_unit_vectors  = get_r_unit_vectors(nside)
local_triad     = get_local_triad(nside,r_unit_vectors)

# Create the magnetic field object
magfield		= MagField(size_box,Npix_box,12345)
# Create the filament population object
#population		= FilPop(nside,Nfil,theta_LH_RMS,size_ratio,size_scale,slope,magfield,1234,dust_template,fixed_distance=False)
population		= FilPop(Nfil,theta_LH_RMS,size_ratio,size_scale,slope,magfield,2345,fixed_distance=False)

TQUmap_total = np.zeros((3,12*nside**2))
for n in range(population.realNfil):
	TQUmap			= FilamentPaint.Paint_Filament(n,nside,local_triad,population,magfield)
	for i in range(3):
		TQUmap_total[i,:] 	+= TQUmap[i,:]

#hp.write_map(output_tqumap,TQUmap_total,nest=False,overwrite=True)