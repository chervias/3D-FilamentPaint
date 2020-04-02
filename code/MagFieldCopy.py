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
		self.Bcube			= np.zeros((self.pixels,self.pixels,self.pixels,3))
		self.Bcube[:,:,:,2] = 1.0
		self.interp_fn		= self.get_interpolator()

def get_interpolator(self):
		real_units				= np.linspace(-0.5*self.size,+0.5*self.size,self.pixels)
		interp_fn			 	= RegularGridInterpolator((real_units,real_units,real_units),self.Bcube,method='linear',fill_value=None)
		return interp_fn