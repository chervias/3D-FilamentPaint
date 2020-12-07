import numpy as np
import healpy as hp

class FilPop:
	def __init__(self,Nfil,theta_LH_RMS,size_ratio,size_scale,slope,Bcube,interpolator,size,seed,alpha,beta,nside,dust_template,beta_template,T_template,ell_limit,fixed_distance=False,fixed_size=False,galactic_plane=False):
		self.Nfil			= Nfil
		#self.realNfil    = Nfil
		self.Bcube		= Bcube
		self.interpolator = interpolator
		self.size = size
		self.fixed_distance	= fixed_distance
		self.fixed_size = fixed_size
		self.nside = nside
		self.galactic_plane = galactic_plane
		self.dust_template = dust_template
		self.beta_template = beta_template
		self.T_template = T_template
		self.max_length		= 1.0
		np.random.seed(seed)
		self.slope			= slope
		self.size_scale		= size_scale
		self.size_ratio		= size_ratio
		self.ell_limit = ell_limit
		self.alpha = alpha
		self.beta = beta
		if theta_LH_RMS == None:
			self.theta_LH_RMS	= None
		elif theta_LH_RMS == -1.0:
			self.theta_LH_RMS	= -1.0
		else:
			self.theta_LH_RMS	= np.radians(theta_LH_RMS)
		self.centers				= self.get_centers()
		self.angles,self.long_vec	= self.get_angles()
		self.sizes					= self.get_sizes()
		self.reject_big_filaments()
		self.fpol0 = self.get_fpol()
		self.beta_array, self.T_array = self.get_beta_T()
		# now delete the objectes according to mask
		self.centers = self.centers[np.logical_not(self.mask),...]
		self.angles = self.angles[np.logical_not(self.mask),...]
		self.long_vec = self.long_vec[np.logical_not(self.mask),...]
		self.sizes = self.sizes[np.logical_not(self.mask),...]
		if self.theta_LH_RMS != -1 and self.theta_LH_RMS != None:
			self.theta_LH = self.theta_LH[np.logical_not(self.mask),...]
			self.psi_LH = self.psi_LH[np.logical_not(self.mask),...]
			self.thetaH = self.thetaH[np.logical_not(self.mask),...]
			self.thetaL = self.thetaL[np.logical_not(self.mask),...]
			self.theta_a = self.theta_a[np.logical_not(self.mask),...]
		self.fpol0 = self.fpol0[np.logical_not(self.mask),...]
		self.beta_array = self.beta_array[np.logical_not(self.mask),...]
		self.T_array = self.T_array[np.logical_not(self.mask),...]
		self.Nfil = self.Nfil - np.sum(self.mask)
	def get_centers(self):
		# recipe to generate random centers
		# Usually this will be done in spherical coordinates
		if self.galactic_plane:
			# Now the centers will be defined by a template map
			# first normalize the map to the total number of filaments
			map_original = hp.read_map(self.dust_template,field=0)
			map_nside = hp.ud_grade(map_original,self.nside)
			# each pixel is mult by C --> Nfil = C*Sum(map)
			C = self.Nfil / np.sum(map_nside)
			number_fil = np.random.poisson(C*map_nside,12*self.nside**2)
			self.realNfil = np.sum(number_fil)
			centers	= np.zeros((self.realNfil,3))
			l_rand	 = np.random.uniform(0.2,1.0,self.realNfil)
			radii_arr = (0.5*self.size)*l_rand**(1.0/3.0)
			counter = 0
			for n in range(12*self.nside**2):
				if number_fil[n] > 0:
					for k in range(number_fil[n]):
						centers[counter,:] = radii_arr[counter]*np.array(hp.pix2vec(self.nside,n,nest=False))
						counter = counter + 1
				else:
					continue
			print('Real number of fil',self.realNfil,'counter',counter)
			# set the number of fil to real fil
			self.Nfil = self.realNfil
			return centers
		elif self.fixed_distance:
			centers	= np.zeros((self.Nfil,3))
			radii_random	= np.ones(self.Nfil) * 0.4*self.size
			phi_random		= 2*np.pi*np.random.uniform(0.0,1.0,self.Nfil)
			theta_random	= np.arccos(1.0 - 2*np.random.uniform(0.0,1.0,self.Nfil))
			centers[:,0]	= radii_random*np.sin(theta_random)*np.cos(phi_random)
			centers[:,1]	= radii_random*np.sin(theta_random)*np.sin(phi_random)
			centers[:,2]	= radii_random*np.cos(theta_random)
			return centers
		else:
			centers	= np.zeros((self.Nfil,3))
			l_rand	 = np.random.uniform(0.2,1.0,self.Nfil)
			u_rand   = np.random.uniform(-1.0,1.0,self.Nfil)
			phi_rand = np.random.uniform(0.0,2.*np.pi,self.Nfil)
			centers[:,0] = (0.5*self.size)*l_rand**(1.0/3.0)*np.sqrt(1. - u_rand**2)*np.cos(phi_rand)
			centers[:,1] = (0.5*self.size)*l_rand**(1.0/3.0)*np.sqrt(1. - u_rand**2)*np.sin(phi_rand)
			centers[:,2] = (0.5*self.size)*l_rand**(1.0/3.0)*u_rand
			return centers
	def get_angles(self):
		angles			= np.zeros((self.Nfil,2))
		# get the euler angles according to the local magnetic field in the center pixel. The hatZ vector of the filament (long axis) follows local B
		local_magfield	= np.array([self.interpolator((self.centers[n,0],self.centers[n,1],self.centers[n,2])) for n in range(self.Nfil)])
		if self.theta_LH_RMS == None:
			# unit vector along the local mag field
			hatZ			= np.array([local_magfield[n,:]/np.linalg.norm(local_magfield[n,:]) for n in range(self.Nfil)])
			# alpha angle
			angles[:,1]		= np.arccos(hatZ[:,2])
			# beta angle
			angles[:,0]		= np.arctan2(hatZ[:,1],hatZ[:,0])
			rhat = np.array([self.centers[n]/np.linalg.norm(self.centers[n]) for n in range(self.Nfil)])
			local_B_proj = np.array([local_magfield[n] - np.dot(local_magfield[n],rhat[n])*rhat[n] for n in range(self.Nfil)])
			filament_vec_proj = np.array([hatZ[n] - np.dot(hatZ[n],rhat[n])*rhat[n] for n in range(self.Nfil)])
			cross = np.array([np.cross(filament_vec_proj[n],local_B_proj[n]) for n in range(self.Nfil)])
			
			return angles,hatZ
		elif self.theta_LH_RMS == -1:
			# we want a unit vector that is ort to center vector
			ort_vec 		= np.array([np.cross(self.centers[n],np.array([1,1,1])) for n in range(self.Nfil)])
			# hatk is the unit vector along the LOS 
			hatk = np.array([self.centers[n]/np.linalg.norm(self.centers[n]) for n in range(self.Nfil)])
			# we rotate the ort vector by a random angle between 0 and 2pi
			phi_angle   = np.random.uniform(0,2*np.pi,self.Nfil)
			ort_vec_rotated		= np.array([ort_vec[n]*np.cos(phi_angle[n]) + np.cross(hatk[n],ort_vec[n])*np.sin(phi_angle[n]) + hatk[n]*np.dot(hatk[n],ort_vec[n])*(1 - np.cos(phi_angle[n])) for n in range(self.Nfil)])
			ort_vec_unit	= np.array([ort_vec_rotated[n]/np.linalg.norm(ort_vec_rotated[n]) for n in range(self.Nfil)])
			# alpha angle
			angles[:,1]		= np.arccos(ort_vec_unit[:,2])
			# beta angle
			angles[:,0]		= np.arctan2(ort_vec_unit[:,1],ort_vec_unit[:,0])
			return angles,ort_vec_unit
		else:
			# unit vector along the local mag field
			hatZ			= np.array([local_magfield[n,:]/np.linalg.norm(local_magfield[n,:]) for n in range(self.Nfil)])
			rhat = np.array([self.centers[n]/np.linalg.norm(self.centers[n]) for n in range(self.Nfil)])
			local_B_proj = np.array([local_magfield[n] - np.dot(local_magfield[n],rhat[n])*rhat[n] for n in range(self.Nfil)])
			# we need a second unit vector hatY perpendicular to hatZ
			random_vectors  = np.random.uniform(0.0,1.0,size=(self.Nfil,3))
			vecY			= np.cross(hatZ,random_vectors)
			hatY			= np.array([vecY[n,:]/np.linalg.norm(vecY[n,:]) for n in range(self.Nfil)])
			# This is in radians
			theta_LH		= np.fabs(np.random.normal(loc=0,scale=self.theta_LH_RMS,size=self.Nfil))
			self.theta_LH = theta_LH
			#theta_LH = np.zeros(self.Nfil)
			phi				= np.random.uniform(0,2*np.pi,self.Nfil)
			# We rotate hatZ around hatY by theta_LH using Rodrigues formula
			hatZprime		= np.array([hatZ[n,:]*np.cos(theta_LH[n]) + np.cross(hatY[n,:],hatZ[n,:])*np.sin(theta_LH[n]) + hatY[n,:]*np.dot(hatY[n,:],hatZ[n,:])*(1 - np.cos(theta_LH[n])) for n in range(self.Nfil)])
			# We rotate hatZprime around hatZ by phi using Rodrigues formula
			hatZprime2		= np.array([hatZprime[n,:]*np.cos(phi[n]) + np.cross(hatZ[n,:],hatZprime[n,:])*np.sin(phi[n]) + hatZ[n,:]*np.dot(hatZ[n,:],hatZprime[n,:])*(1 - np.cos(phi[n])) for n in range(self.Nfil)])
			filament_vec_proj = np.array([hatZprime2[n] - np.dot(hatZprime2[n],rhat[n])*rhat[n] for n in range(self.Nfil)])
			# I calculate the psi_LH angles, with cw/ccw sign
			self.psi_LH = np.array([np.arctan2(np.dot(rhat[n],np.cross(filament_vec_proj[n],local_B_proj[n])),np.dot(filament_vec_proj[n],local_B_proj[n])) for n in range(self.Nfil)])
			self.thetaH = np.array([np.arccos(np.dot(local_magfield[n],rhat[n])/np.linalg.norm(local_magfield[n])) for n in range(self.Nfil)])
			self.thetaL = np.array([np.arccos(np.dot(hatZprime2[n],rhat[n])/np.linalg.norm(hatZprime2[n])) for n in range(self.Nfil)])
			# Now hatZprime2 is the direction of the long axis of the filament
			norm_hatZprime2	= np.linalg.norm(hatZprime2,axis=1)
			# beta angle
			angles[:,1]		= np.arccos(hatZprime2[:,2]/norm_hatZprime2)
			# alpha angle
			angles[:,0]		= np.arctan2(hatZprime2[:,1],hatZprime2[:,0])
			return angles,hatZprime2
	def get_sizes(self):
		# The sizes will be the ellipsoid semi axes a,b,c with a=b<c
		sizes			= np.zeros((self.Nfil,3))
		if self.fixed_size:
			c_semiaxis = self.size_scale*np.ones(self.Nfil)
		else:
			c_semiaxis		= self.size_scale*(1.0+np.random.pareto(self.slope-1,size=self.Nfil))
		#a_semiaxis		= self.size_scale*(1.0+np.random.pareto(self.slope-1,size=self.Nfil))
		#sizes[:,2]		= np.clip(c_semiaxis,0,self.max_length)
		sizes[:,2]		= c_semiaxis
		sizes[:,0]		= self.size_ratio*sizes[:,2]
		sizes[:,1]		= self.size_ratio*sizes[:,2]
		return sizes
	def get_fpol(self):
		# create a beta distribution for fpol0 
		fpol0 = np.random.beta(self.alpha,self.beta,size=self.Nfil)
		return fpol0
	def get_beta_T(self):
		beta_map_original = hp.read_map(self.beta_template,field=(0,1))
		beta_map_nside = hp.ud_grade(beta_map_original,self.nside)
		# load correlation
		correlation = np.load('/home/chervias/CMB-work/Filaments/3dfilament-project-healpix/output/ThermalDustCorrelation.npz')
		beta_array = np.zeros((self.realNfil))
		#T_array = np.zeros((self.realNfil))
		for n in range(self.realNfil):
			# determine which pixel corresponds to center
			pixel = hp.vec2pix(self.nside,self.centers[n,0],self.centers[n,1],self.centers[n,2])
			beta_array[n] = np.random.normal(loc=beta_map_nside[0][pixel],scale=0.9*beta_map_nside[1][pixel],size=1)
		# correlated Tdust
		T_array = np.interp(beta_array,correlation['beta_array'], correlation['fitted_T'])
		return beta_array,T_array
	def reject_big_filaments(self):
		self.theta_a = np.array([np.sqrt((self.sizes[n,2]*np.cos(np.arctan(np.tan(self.thetaL[n])/self.size_ratio)))**2 + (self.sizes[n,1]*np.sin(np.arctan(np.tan(self.thetaL[n])/self.size_ratio)))**2)/np.linalg.norm(self.centers[n]) for n in range(self.Nfil)])
		self.mask =  np.pi / self.theta_a < self.ell_limit

#---------------------------------------------------------------------