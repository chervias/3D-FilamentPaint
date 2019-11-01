import numpy as np
from scipy.interpolate import RegularGridInterpolator

class MagField_noB:
	def __init__(self,size,filename):
		# size is a physical size
		self.size			= size
		self.Bcube			= self.get_MagField()
		self.interp_fn		= self.get_interpolator()
	
	def get_MagField(self):
		Bcube					= np.load(self.filename,)
		return Bcube
	def get_interpolator(self):
		real_units				= np.linspace(-0.5*self.size,+0.5*self.size,self.pixels)
		interp_fn			 	= RegularGridInterpolator((real_units,real_units,real_units),self.Bcube,method='linear',fill_value=None)
		return interp_fn
