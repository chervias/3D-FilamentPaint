import FilamentPaint

import sys
sys.path.insert(0,'../src')
from Sky import Sky
from MagField import MagField
from FilPop import FilPop
from Filament import Filament

nside			= 512
Npix_box		= 128
Nfil			= 1
theta_LH_RMS	= 10. # in degrees
size_box		= 1500.0 # physical size box
size_scale		= 0.7
size_ratio		= 0.25

# Create the sky object
sky				= Sky(nside)
# Create the magnetic field object
magfield		= MagField(size_box,Npix_box,12345)
# Create the filament population object
population		= FilPop(Nfil,theta_LH_RMS,size_ratio,size_scale,magfield,fixed_distance=True)

sizes 			= population.sizes[0]
angles 			= population.angles[0]
center 			= population.centers[0]
filament		= Filament(center,sizes,angles)

result	= FilamentPaint.Paint_Filament(0,sky,population,magfield)
print(result)
print('------------------------')
print(filament.xyz_edges_vectors)
