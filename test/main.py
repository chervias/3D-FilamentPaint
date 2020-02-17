import FilamentPaint
import numpy as np
import healpy as hp
from mpi4py import MPI
import sys
sys.path.insert(0,'code')
from Sky import Sky
from MagField import MagField
from FilPop import FilPop

output_tqumap	= 'test.fits'
nside			= 1024
Npix_box		= 256
theta_LH_RMS	= -1 # in degrees
size_scale		= 0.1
size_ratio		= 0.25
Nfil			= 100
size_box		= 400.0 # physical size box
slope			= 2.4

# Create the sky object
sky				= Sky(nside)
# Create the magnetic field object
magfield		= MagField(size_box,Npix_box,12345)
# Create the filament population object
population		= FilPop(nside,Nfil,theta_LH_RMS,size_ratio,size_scale,slope,magfield,1234,fixed_distance=True)

for n in range(Nfil):
	a = FilamentPaint.Paint_Filament(n,sky,population,magfield)
