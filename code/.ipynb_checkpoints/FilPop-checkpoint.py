import numpy as np

class FilPop:
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
		self.centers				= self.get_centers()	
		self.angles,self.long_vec	= self.get_angles()
		self.sizes					= self.get_sizes()
	def get_centers(self):
		centers	= np.zeros((self.Nfil,3))
		# recipe to generate random centers
		# Usually this will be done in spherical coordinates
		if self.fixed_distance:
			radii_random	= np.ones(self.Nfil) * 0.4*self.magfield.size
			phi_random		= 2*np.pi*np.random.uniform(0.0,1.0,self.Nfil)
			theta_random	= np.arccos(1.0 - 2*np.random.uniform(0.0,1.0,self.Nfil))
			centers[:,0]	= radii_random*np.sin(theta_random)*np.cos(phi_random)
			centers[:,1]	= radii_random*np.sin(theta_random)*np.sin(phi_random)
			centers[:,2]	= radii_random*np.cos(theta_random)
		else:
			l_rand	 = np.random.uniform(0.0,1.0,self.Nfil)
			u_rand   = np.random.uniform(-1.0,1.0,self.Nfil)
			phi_rand = np.random.uniform(0.0,2.*np.pi,self.Nfil)
			centers[:,0] = (0.5*self.magfield.size)*l_rand**(1.0/3.0)*np.sqrt(1. - u_rand**2)*np.cos(phi_rand)
			centers[:,1] = (0.5*self.magfield.size)*l_rand**(1.0/3.0)*np.sqrt(1. - u_rand**2)*np.sin(phi_rand)
			centers[:,2] = (0.5*self.magfield.size)*l_rand**(1.0/3.0)*u_rand
		return centers
	def get_angles(self):
		angles			= np.zeros((self.Nfil,2))
		# get the euler angles according to the local magnetic field in the center pixel. The hatZ vector of the filament (long axis) follows local B
		local_magfield	= np.array([self.magfield.interp_fn((self.centers[n,0],self.centers[n,1],self.centers[n,2])) for n in range(self.Nfil)])
		if self.theta_LH_RMS == None:
			# unit vector along the local mag field
			hatZ			= np.array([local_magfield[n,:]/np.linalg.norm(local_magfield[n,:]) for n in range(self.Nfil)])
			# alpha angle
			angles[:,1]		= np.arccos(hatZ[:,2])
			# beta angle
			angles[:,0]		= np.arctan2(hatZ[:,1],hatZ[:,0])
			return angles,hatZ
		elif self.theta_LH_RMS == -1:
			# we want a unit vector that is ort to center vector
			ort_vec 		= np.array([np.cross(self.centers[n],np.array([1,1,1])) for n in range(self.Nfil)])
			ort_vec_unit	= np.array([ort_vec[n,:]/np.linalg.norm(ort_vec[n,:]) for n in range(self.Nfil)])
			# alpha angle
			angles[:,1]		= np.arccos(ort_vec_unit[:,2])
			# beta angle
			angles[:,0]		= np.arctan2(ort_vec_unit[:,1],ort_vec_unit[:,0])
			return angles,ort_vec_unit
		else:
			# unit vector along the local mag field
			hatZ			= np.array([local_magfield[n,:]/np.linalg.norm(local_magfield[n,:]) for n in range(self.Nfil)])
			# we need a second unit vector hatY perpendicular to hatZ
			random_vectors  = np.random.uniform(0.0,1.0,shape=())
			vecY			= np.cross(hatZ,np.array([0,1,0]))
			hatY			= np.array([vecY[n,:]/np.linalg.norm(vecY[n,:]) for n in range(self.Nfil)])
			# This is in radians
			theta_LH		= np.fabs(np.random.normal(0,self.theta_LH_RMS,self.Nfil))
			#theta_LH = np.zeros(self.Nfil)
			phi				= np.random.uniform(0,2*np.pi,self.Nfil)
			#phi = np.zeros(self.Nfil)
			# We rotate hatZ around hatY by theta_LH using Rodrigues formula
			hatZprime		= np.array([hatZ[n,:]*np.cos(theta_LH[n]) + np.cross(hatY[n,:],hatZ[n,:])*np.sin(theta_LH[n]) + hatY[n,:]*np.dot(hatY[n,:],hatZ[n,:])*(1 - np.cos(theta_LH[n])) for n in range(self.Nfil)])
			# We rotate hatZprime around hatZ by phi using Rodrigues formula
			hatZprime2		= np.array([hatZprime[n,:]*np.cos(phi[n]) + np.cross(hatZ[n,:],hatZprime[n,:])*np.sin(phi[n]) + hatZ[n,:]*np.dot(hatZ[n,:],hatZprime[n,:])*(1 - np.cos(phi[n])) for n in range(self.Nfil)])
			# Now hatZprime2 is the direction of the long axis of the filament
			norm_hatZprime2	= np.linalg.norm(hatZprime2,axis=1)
			# alpha angle
			angles[:,1]		= np.arccos(hatZprime2[:,2]/norm_hatZprime2)
			# beta angle
			angles[:,0]		= np.arctan2(hatZprime2[:,1],hatZprime2[:,0])
			return angles,hatZprime2
	def get_sizes(self):
		# The sizes will be the ellipsoid semi axes a,b,c with a=b<c
		sizes			= np.zeros((self.Nfil,3))
		c_semiaxis		= self.size_scale*(1.0+np.random.pareto(self.slope-1,size=self.Nfil))
		#sizes[:,2]		= np.clip(c_semiaxis,0,self.max_length)
		sizes[:,2]		= c_semiaxis
		sizes[:,0]		= self.size_ratio*sizes[:,2]
		sizes[:,1]		= self.size_ratio*sizes[:,2]
		return sizes

#---------------------------------------------------------------------