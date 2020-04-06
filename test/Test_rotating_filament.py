import FilamentPaint
import numpy as np
import healpy as hp
import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
from mpi4py import MPI
import sys
sys.path.insert(0,'../code')
from Sky import *
from MagField import MagField
from FilPop import FilPop
import matplotlib.animation as animation

nside			= 512
Npix_box		= 256
theta_LH_RMS	= -1 # in degrees, if -1 the filaments are perpendicular to the LOS
size_scale		= 0.025 # size scale
size_ratio		= 0.25
slope			= 2.4
Nfil			= 1
size_box		= 400.0

Nangles = 25

class FilPopEmpty:
	def __init__(self,Nfil,theta_LH_RMS,size_ratio,size_scale,slope,magfield,seed,fixed_distance=False):
		self.Nfil			= Nfil
		self.magfield		= magfield
		self.fixed_distance	= fixed_distance
		self.max_length		= 1.0
		np.random.seed(seed)
		self.slope			= slope
		self.size_scale		= size_scale
		self.size_ratio		= size_ratio
		if theta_LH_RMS == None:
			self.theta_LH_RMS	= None
		elif theta_LH_RMS == -1.0:
			self.theta_LH_RMS	= -1.0
		else:
			self.theta_LH_RMS	= np.radians(theta_LH_RMS)
		self.centers				= np.zeros((Nfil,3))
		self.angles,self.long_vec	= np.zeros((Nfil,2)),np.zeros((Nfil,3))
		self.sizes					= np.zeros((Nfil,3))
		
r_unit_vectors  = get_r_unit_vectors(nside)
local_triad     = get_local_triad(nside,r_unit_vectors)

# Modify the magfield
angle_arr = np.linspace(0,2*np.pi,Nangles)

center = np.array([75,75,0])
pix_id = hp.vec2pix(nside,center[0],center[1],center[2],nest=False)
XYZ = local_triad[pix_id]
magfield		= MagField(400.0,256,23456)

pop_arr = []
khat = np.array([0.435,-0.235,1.224])/np.linalg.norm(np.array([0.435,-0.235,1.224]))
ort_vec = np.cross(khat,np.array([1,2,2]))

# modify the pop object
for angle in angle_arr:
	population = FilPopEmpty(Nfil,theta_LH_RMS,size_ratio,size_scale,slope,magfield,1234,fixed_distance=True)
	population.centers[0] = center
	ort_vec_rot = ort_vec*np.cos(angle) + np.cross(khat,ort_vec)*np.sin(angle) + khat * np.dot(khat,ort_vec) *(1.0 - np.cos(angle))
	population.angles[0,1]= np.arccos(ort_vec_rot[2])
	population.angles[0,0]= np.arctan2(ort_vec_rot[1],ort_vec_rot[0])
	population.sizes[0] = np.array([2,2,10])
	pop_arr.append(population)

tqu_arr = []
for population in pop_arr:
    tqu_total = FilamentPaint.Paint_Filament(0,nside,local_triad,population,magfield)
    tqu_arr.append(tqu_total)
#(h,w) = pl.figaspect(1.0)
#fig = pl.figure(figsize=(h,w))
fig = pl.figure()

ims = []
for n in range(Nangles):
	tqu_map = tqu_arr[n]
	#almT,almE,almB = hp.map2alm(tqu_map,pol=True)
	#E_map = hp.alm2map(almE,nside)
	#B_map = hp.alm2map(almB,nside)
	T = hp.gnomview(tqu_map[0],rot=(45,0),reso=10.0,nest=False,notext=True,title='',cbar=False,return_projected_map=True,no_plot=True,min=-1,max=1)
	ax_T = pl.gca() ; ax_T.set_aspect('equal');ax_T.set_xticklabels([]);ax_T.set_yticklabels([])
	im = pl.imshow(T,vmin=-1,vmax=1,animated=True)
	ims.append([im])

ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=1000)
ani.save('dynamic_images.mp4')