import sys
sys.path.insert(0,'src')
import numpy as np
import healpy as hp
import matplotlib.pyplot as pl
from Functions import *
from Sky import Sky
from MagField_noB import MagField_noB
from FilPop import FilPop
import matplotlib.pyplot as pl

output_tqumap	= 'test_ns512_Nf10k_maa50_sr0p25_sl4p6.fits'
nside			= 512
Npix_box		= 256
theta_LH_RMS	= 50.0 # in degrees
size_scale		= 1.0
size_ratio		= 0.25
Nfil			= 10000
size_box		= 1500.0 # physical size box
# Create the sky object
sky				= Sky(nside)
# Create the magnetic field object
magfield		= MagField_noB(size_box,'magfield_Bcube.npy')
# Create the filament population object
population		= FilPop(Nfil,theta_LH_RMS,size_ratio,size_scale,magfield)

angsize			= np.zeros((Nfil,2))

for n in range(Nfil):
	#which pixel corresponds to the center of the filament
	#print(population.centers[n].shape)
	ipix 			= hp.vec2pix(nside,population.centers[n,0],population.centers[n,1],population.centers[n,2],nest=False)
	LOS_angle		= np.arccos(np.dot(population.long_vec[n],sky.r_unit_vectors[n]))
	Radial_distance	= np.linalg.norm(population.centers[n])
	angsize[n,0]	= np.sqrt(np.sin(LOS_angle)**2 + size_ratio**2*np.cos(LOS_angle)**2) * 2*population.sizes[n,2] / Radial_distance

pl.hist(2*np.pi/angsize[:,0],np.linspace(2,8*nside,200))
pl.gca().axvline(3*nside,c='red')
pl.gca().set_xlabel(r'multipole $\ell$')
pl.gca().set_xlim(0,8*nside)


#pl.show()
pl.savefig('histogram_ang_sizes.pdf',format='pdf')
