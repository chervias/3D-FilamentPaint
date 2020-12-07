import Bpowspec
import numpy as np
from scipy.interpolate import RegularGridInterpolator

def get_MagField(size,pixels,seed,method,path):
	# WARNING !!!!!!!
	# Kevin's code will output a B cube with indices Bcube[iz,iy,ix], so we have to transpose to put it in the format
	# [ix,iy,iz]
	# !!!!!!!!!!!!!!!
	# Load the large scale Bcube
	
	if method==1:
		Bcube_ls = np.load(path)['Bcube'] # this is in [iz,iy,ix] order
		size3d 					= np.array([size,size,size])
		N 						= np.array([pixels,pixels,pixels],dtype=np.int32)
		Delta 					= size3d/N
		Deltak 					= Bpowspec.Delta2k(Delta,N)
		kmax 					= np.amax(N*Deltak)
		k 						= np.linspace(0,kmax, 1000);
		# I create a P(k)
		#filter = 0.5+0.5*np.tanh(0.015*(k - 900))
		#Pk_small = 1e0*k**-4*filter
		#Pk_small[0] = 0.0
		# This is for the cube with 20kpc per side
		filter =  0.5+0.5*np.tanh(1.0*(k - 15))
		Pk_small = 2e-3*k**-4*filter
		Pk_small[0] = 0.0
		Bharmx,Bharmy,Bharmz 	= Bpowspec.Pk2harm(k,Pk_small,N,kmax,Deltak,seed)
		Bcube					= np.zeros((pixels,pixels,pixels,3))
		Bcube[:,:,:,0] 			= np.transpose(Bpowspec.harm2map(Bharmx,Delta),axes=(2,1,0)) + np.transpose(Bcube_ls[:,:,:,0],axes=(2,1,0))
		Bcube[:,:,:,1] 			= np.transpose(Bpowspec.harm2map(Bharmy,Delta),axes=(2,1,0)) + np.transpose(Bcube_ls[:,:,:,1],axes=(2,1,0))
		Bcube[:,:,:,2] 			= np.transpose(Bpowspec.harm2map(Bharmz,Delta),axes=(2,1,0)) + np.transpose(Bcube_ls[:,:,:,2],axes=(2,1,0))
		return Bcube
	if method==2:
		# WARNING !!!!!!!
		# Kevin's code will output a B cube with indices Bcube[iz,iy,ix], so we have to transpose to put it in the format
		# [ix,iy,iz]
		# !!!!!!!!!!!!!!!
		size3d 					= np.array([size,size,size])
		N 						= np.array([pixels,pixels,pixels],dtype=np.int32)
		Delta 					= size3d/N
		Deltak 					= Bpowspec.Delta2k(Delta,N)
		kmax 					= np.amax(N*Deltak)
		k 						= np.linspace(0,kmax, 500);
		#Pk 						= np.exp(-k**2/2/(kmax/40)**2)
		Pk 						= np.exp(-k**2/2/(kmax/100)**2)
		Bharmx,Bharmy,Bharmz 	= Bpowspec.Pk2harm(k,Pk,N,kmax,Deltak,seed)
		Bcube					= np.zeros((pixels,pixels,pixels,3))
		Bcube[:,:,:,0] 			= np.transpose(Bpowspec.harm2map(Bharmx,Delta),axes=(2,1,0))
		Bcube[:,:,:,1] 			= np.transpose(Bpowspec.harm2map(Bharmy,Delta),axes=(2,1,0))
		Bcube[:,:,:,2] 			= np.transpose(Bpowspec.harm2map(Bharmz,Delta),axes=(2,1,0))
		return Bcube
	if method==3:
		# WARNING !!!!!!!
		# Kevin's code will output a B cube with indices Bcube[iz,iy,ix], so we have to transpose to put it in the format
		# [ix,iy,iz]
		# !!!!!!!!!!!!!!!
		size3d 					= np.array([size,size,size])
		N 						= np.array([pixels,pixels,pixels],dtype=np.int32)
		Delta 					= size3d/N
		Deltak 					= Bpowspec.Delta2k(Delta,N)
		kmax 					= np.amax(N*Deltak)
		k 						= np.linspace(0,kmax, 1000);
		import healpy as hp
		beam = hp.gauss_beam(np.radians(0.35),lmax=999,pol=False)
		Pk = k**-4*beam**2
		Pk[0] = 0.0
		Pk[-1] = 0.0
		Bharmx,Bharmy,Bharmz 	= Bpowspec.Pk2harm(k,Pk,N,kmax,Deltak,seed)
		Bcube					= np.zeros((pixels,pixels,pixels,3))
		Bcube[:,:,:,0] 			= np.transpose(Bpowspec.harm2map(Bharmx,Delta),axes=(2,1,0))
		Bcube[:,:,:,1] 			= np.transpose(Bpowspec.harm2map(Bharmy,Delta),axes=(2,1,0))
		Bcube[:,:,:,2] 			= np.transpose(Bpowspec.harm2map(Bharmz,Delta),axes=(2,1,0))
		return Bcube

def get_interpolator(size,pixels,Bcube):
	real_units				= np.linspace(-0.5*size,+0.5*size,pixels)
	interp_fn			 	= RegularGridInterpolator((real_units,real_units,real_units),Bcube,method='linear',fill_value=None)
	return interp_fn