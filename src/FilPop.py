import numpy as np

class FilPop:
	def __init__(self,Nfil,theta_LH_RMS,magfield):
		self.Nfil			= Nfil
		self.magfield		= magfield
		self.max_length		= 10.0
		self.size_ratio		= 0.25
		if theta_LH_RMS == None:
			self.theta_LH_RMS	= None
		else:
			self.theta_LH_RMS	= np.radians(theta_LH_RMS)
		self.centers		= self.get_centers()	
		self.angles			= self.get_angles()
		self.sizes			= self.get_sizes()
	def get_centers(self):
		centers	= np.zeros((self.Nfil,3))
		# recipe to generate random centers
		# Usually this will be done in spherical coordinates
		# the center cannot leave the box within where B is defined
		# maximum radial distance of a center is 0.5*size - 5*max_length
		radii_random	= np.random.uniform((0.05*self.magfield.size)**3,(0.5*self.magfield.size - 5*self.max_length)**3,self.Nfil)**(1./3.)
		phi_random		= 2*np.pi*np.random.uniform(0.0,1.0,self.Nfil)
		theta_random	= np.arccos(1.0 - 2*np.random.uniform(0.0,1.0,self.Nfil))
		centers[:,0]	= radii_random*np.sin(theta_random)*np.cos(phi_random)
		centers[:,1]	= radii_random*np.sin(theta_random)*np.sin(phi_random)
		centers[:,2]	= radii_random*np.cos(theta_random)
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
			return angles
		else:
			# unit vector along the local mag field
			hatZ			= np.array([local_magfield[n,:]/np.linalg.norm(local_magfield[n,:]) for n in range(self.Nfil)])
			# we need a second unit vector hatY perpendicular to hatZ
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
			return angles
	def get_sizes(self):
		# The sizes will be the ellipsoid semi axes a,b,c with a=b<c
		sizes			= np.zeros((self.Nfil,3))
		c_semiaxis		= 1 + np.random.pareto(4.6-1,size=self.Nfil)
		sizes[:,2]		= np.clip(c_semiaxis,0,self.max_length)
		sizes[:,0]		= self.size_ratio*sizes[:,2]
		sizes[:,1]		= self.size_ratio*sizes[:,2]
		return sizes

class FilPop_1fil:
	# This is to test with a single filament
	def __init__(self,center,sizes,theta_LH,magfield):
		self.magfield	= magfield
		self.theta_LH	= np.array(theta_LH)
		self.centers	= np.array([center])
		self.angles		= self.get_angles()
		self.sizes		= np.array([sizes])
	def get_angles(self):
		angles			= np.zeros((1,2))
		# get the euler angles according to the local magnetic field in the center pixel. The hatZ vector of the filament (long axis) follows local B
		local_magfield	= self.magfield.interp_fn((self.centers[0,0],self.centers[0,1],self.centers[0,2]))
		# unit vector along the local mag field
		hatZ			= local_magfield/np.linalg.norm(local_magfield)
		hatX			= np.array([1,0,0])
		#print(hatZ)
		# This is in radians
		theta_LH		= np.radians(self.theta_LH)
		phi				= 0.0
		# We rotate hatZ around hatY by theta_LH using Rodrigues formula
		hatZprime		= hatZ*np.cos(theta_LH) + np.cross(hatX,hatZ)*np.sin(theta_LH) + hatX*np.dot(hatX,hatZ)*(1 - np.cos(theta_LH))
		# We rotate hatZprime around hatZ by phi using Rodrigues formula
		hatZprime2		= hatZprime
		# Now hatZprime2 is the direction of the long axis of the filament
		norm_hatZprime2	= np.linalg.norm(hatZprime2)
		# alpha angle
		angles[0,1]		= np.arccos(hatZprime2[2]/norm_hatZprime2)
		# beta angle
		angles[0,0]		= np.arctan2(hatZprime2[1],hatZprime2[0])
		return angles
