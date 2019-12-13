import FilamentPaint

import healpy as hp
import matplotlib.pyplot as pl

import sys
sys.path.insert(0,'../code')
from Sky import Sky
from MagField import MagField
from FilPop import FilPop
from Filament import Filament
from Functions import *

nside			= 512
Npix_box		= 128
Nfil			= 1
theta_LH_RMS	= None # in degrees
size_box		= 1800.0 # physical size box
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

pix_filament	= filament.do_query_polygon(sky)
N_pix_filament	= len(pix_filament)
for n in range(N_pix_filament):
	idx_pixel			= pix_filament[n]
	r_distances			= distances(filament,sky,idx_pixel)
	# check that they are 2
	if len(r_distances) != 2:
		print('Error: something went wrong in the intersection of rays and faces')
		exit()
	r1	= r_distances[0]
	r2	= r_distances[1]		
	int_t,int_q,int_u		= integrator(filament,sky,idx_pixel,magfield,r1,r2,5)
	print("ipix=%i (%.5E,%.5E,%.5E) "%(idx_pixel,int_t,int_q,int_u))

print('---------------------------------')

TQUmap			= FilamentPaint.Paint_Filament(0,sky,population,magfield)
#tqumap			= np.array(TQUmap)

for n in range(N_pix_filament):
	idx_pixel			= pix_filament[n]
	print("ipix=%i (%.5E,%.5E,%.5E) "%(idx_pixel,TQUmap[0,idx_pixel],TQUmap[1,idx_pixel],TQUmap[2,idx_pixel]))

