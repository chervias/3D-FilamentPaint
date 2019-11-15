import sys
sys.path.insert(0,'src')
import numpy as np
import healpy as hp
import matplotlib.pyplot as pl
from Functions import *
from Sky import Sky
from MagField import MagField
from FilPop import FilPop
import matplotlib.pyplot as pl

nside			= 512
Npix_box		= 256
theta_LH_RMS	= None # in degrees
size_ratio		= 0.25
Nfil			= 100000
size_box		= 1500.0 # physical size box
# Create the sky object
sky				= Sky(nside)
# Create the magnetic field object
magfield		= MagField(size_box,Npix_box,12345)
# Create the filament population object
population		= FilPop(Nfil,theta_LH_RMS,size_ratio,1,magfield,fixed_distance=True)
population2		= FilPop(Nfil,theta_LH_RMS,size_ratio,0.7,magfield,fixed_distance=True)

ax = pl.gca()

angsize		= np.zeros((Nfil,2))
for n in range(Nfil):
	#which pixel corresponds to the center of the filament
	#print(population.centers[n].shape)
	ipix 			= hp.vec2pix(nside,population.centers[n,0],population.centers[n,1],population.centers[n,2],nest=False)
	LOS_angle		= np.arccos(np.dot(population.long_vec[n],sky.r_unit_vectors[n]))
	Radial_distance	= np.linalg.norm(population.centers[n])
	angsize[n,0]	= np.sqrt(np.sin(LOS_angle)**2 + size_ratio**2*np.cos(LOS_angle)**2) * 2*population.sizes[n,2] / Radial_distance

#ax.hist(2*np.pi/angsize[:,0],np.linspace(2,8*nside,200),histtype='step',label='Size scale = 1')
ax.hist(angsize[:,0],np.linspace(0,0.02,200),histtype='step',label='Size scale = 1')
	
angsize		= np.zeros((Nfil,2))
for n in range(Nfil):
	#which pixel corresponds to the center of the filament
	#print(population.centers[n].shape)
	ipix 			= hp.vec2pix(nside,population2.centers[n,0],population2.centers[n,1],population2.centers[n,2],nest=False)
	LOS_angle		= np.arccos(np.dot(population2.long_vec[n],sky.r_unit_vectors[n]))
	Radial_distance	= np.linalg.norm(population2.centers[n])
	angsize[n,0]	= np.sqrt(np.sin(LOS_angle)**2 + size_ratio**2*np.cos(LOS_angle)**2) * 2*population2.sizes[n,2] / Radial_distance

#ax.hist(2*np.pi/angsize[:,0],np.linspace(2,8*nside,200),histtype='step',color='green',label='Size scale = 0.7')
ax.hist(angsize[:,0],np.linspace(0,0.02,200),histtype='step',color='green',label='Size scale = 0.7')

ax.set_title('%i filaments'%Nfil)
#ax.axvline(3*nside,c='red')
ax.set_xlabel(r'Angular size [rad]')
#ax.set_xlim(0,8*nside)
ax.legend(loc='best')

#pl.show()
pl.savefig('histogram_ang_sizes.pdf',format='pdf')
