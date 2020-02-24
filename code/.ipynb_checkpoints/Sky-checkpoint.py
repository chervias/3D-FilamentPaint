import numpy as np
import healpy as hp

class Sky:
	def __init__(self,nside):
		self.nside								= nside
		self.Tmap								= np.zeros(12*self.nside**2)
		self.Qmap								= np.zeros(12*self.nside**2)
		self.Umap								= np.zeros(12*self.nside**2)
		self.mask								= np.zeros(12*self.nside**2)
		self.r_unit_vectors						= self.get_r_unit_vectors()
		self.local_triad						= self.get_local_triad()
			
	def get_r_unit_vectors(self):
		#get the r unit vector for each pixel in a ring pixelization and the angles
		Npix	= 12*self.nside**2
		r_unit_vectors			= np.zeros((Npix,3))
		result 					= hp.pix2vec(self.nside,range(Npix),nest=False)
		r_unit_vectors[:,0]		= np.array(result[0])
		r_unit_vectors[:,1]		= np.array(result[1])
		r_unit_vectors[:,2]		= np.array(result[2])
		return r_unit_vectors
	
	def get_local_triad(self):
		# Return local coordinate system at specified vector nhat, according to Healpix polarization convention
		# note: fails for nhat is zhat0 = North Pole
		Npix				= 12*self.nside**2
		# the indices are (N_pix,Ntriad_vectors,Ndimensions)
		local_triad			= np.zeros((Npix,3,3))
		nhat				= self.r_unit_vectors
		zhat0 				= np.array([0,0,1])
		local_triad[:,2,:] 	= nhat
		tmp 				= -zhat0 + np.einsum('i,ij->ij',np.dot(nhat,zhat0),nhat)
		local_triad[:,0,:] 	= tmp / np.linalg.norm(tmp,axis=1)[:,None]
		local_triad[:,1,:] 	= np.cross(local_triad[:,2],local_triad[:,0])
		return local_triad
		
	def save_sky(self,name):
		hp.write_map(name,[self.Tmap,self.Qmap,self.Umap],nest=False,overwrite=True)

class SkyAux:
	def __init__(self,nside):
		self.nside								= nside
		self.Tmap								= None
		self.Qmap								= None
		self.Umap								= None
		self.mask								= None
		self.r_unit_vectors						= None
		self.local_triad						= None