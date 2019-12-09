import Bpowspec
import numpy as np
from scipy.interpolate import RegularGridInterpolator

class MagField:
	def __init__(self,size,pixels,seed):
		# size is a physical size
		self.size			= size
		# pixels is the number of pixels on each dimension
		self.pixels			= pixels
		self.seed			= seed
		self.Bcube			= self.get_MagField()
		self.interp_fn		= self.get_interpolator()
	
	def get_MagField(self):
		size3d 					= np.array([self.size,self.size,self.size])
		N 						= np.array([self.pixels,self.pixels,self.pixels],dtype=np.int32)
		Delta 					= size3d/N
		Deltak 					= Bpowspec.Delta2k(Delta,N)
		kmax 					= np.amax(N*Deltak)
		k 						= np.linspace(0,kmax, 500);
		Pk 						= np.exp(-k**2/2/(kmax/40)**2)
		kx,ky,kz 				= Bpowspec.kvecs(N,Deltak)
		Bharmx,Bharmy,Bharmz 	= Bpowspec.Pk2harm(k,Pk,N,kmax,Deltak,self.seed)
		Bcube					= np.zeros((self.pixels,self.pixels,self.pixels,3))
		Bcube[:,:,:,0] 			= Bpowspec.harm2map(Bharmx,Delta)
		Bcube[:,:,:,1] 			= Bpowspec.harm2map(Bharmy,Delta)
		Bcube[:,:,:,2] 			= Bpowspec.harm2map(Bharmz,Delta)
		return Bcube
	
	def get_interpolator(self):
		real_units				= np.linspace(-0.5*self.size,+0.5*self.size,self.pixels)
		interp_fn			 	= RegularGridInterpolator((real_units,real_units,real_units),self.Bcube,method='linear',fill_value=None)
		return interp_fn

class MagField_1fil:
	# This is to test with a single filament
	def __init__(self,size,pixels,seed,direction):
		# size is a physical size
		self.size			= size
		# pixels is the number of pixels on each dimension
		self.pixels			= pixels
		self.direction		= direction
		self.seed			= seed
		self.Bcube			= self.get_MagField()
		self.interp_fn		= self.get_interpolator()
	
	def get_MagField(self):
		Bcube					= np.zeros((self.pixels,self.pixels,self.pixels,3))
		
		if self.direction=='+z':
			Bcube[:,:,:,0]			= np.zeros((self.pixels,self.pixels,self.pixels))
			Bcube[:,:,:,1]			= np.zeros((self.pixels,self.pixels,self.pixels))
			Bcube[:,:,:,2]			= np.ones((self.pixels,self.pixels,self.pixels))
		if self.direction=='-z':
			Bcube[:,:,:,0]			= np.zeros((self.pixels,self.pixels,self.pixels))
			Bcube[:,:,:,1]			= np.zeros((self.pixels,self.pixels,self.pixels))
			Bcube[:,:,:,2]			= -1.0*np.ones((self.pixels,self.pixels,self.pixels))
		if self.direction=='+y':
			Bcube[:,:,:,0]			= np.zeros((self.pixels,self.pixels,self.pixels))
			Bcube[:,:,:,1]			= np.ones((self.pixels,self.pixels,self.pixels))
			Bcube[:,:,:,2]			= np.zeros((self.pixels,self.pixels,self.pixels))
		if self.direction=='-y':
			Bcube[:,:,:,0]			= np.zeros((self.pixels,self.pixels,self.pixels))
			Bcube[:,:,:,1]			= -1*np.ones((self.pixels,self.pixels,self.pixels))
			Bcube[:,:,:,2]			= np.zeros((self.pixels,self.pixels,self.pixels))
		if self.direction=='45deg':
			Bcube[:,:,:,0]			= np.zeros((self.pixels,self.pixels,self.pixels))
			Bcube[:,:,:,1]			= np.ones((self.pixels,self.pixels,self.pixels))
			Bcube[:,:,:,2]			= -1*np.ones((self.pixels,self.pixels,self.pixels))
		return Bcube
	
	def get_interpolator(self):
		real_units				= np.linspace(-0.5*self.size,+0.5*self.size,self.pixels)
		interp_fn			 	= RegularGridInterpolator((real_units,real_units,real_units),self.Bcube,method='linear',fill_value=None)
		return interp_fn
