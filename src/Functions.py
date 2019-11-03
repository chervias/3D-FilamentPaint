import numpy as np
import healpy as hp
from scipy.special.orthogonal import p_roots
from Filament import Filament

def distances(filament,sky,idx_pixel):
	#calculate the intersection distance between the r*hat(r) vector and the plane defined by each of the faces
	radii_intersect	= np.zeros(6)
	id_faces		= []
	for k in range(6):
		# radius of intersection between vector r*hat(r) and the plane of face k
		bottom	= np.dot(filament.xyz_normal[k],sky.r_unit_vectors[idx_pixel])
		if bottom == 0.0:
			radii_intersect[k]		= np.nan
		else:
			radii_intersect[k]			= np.dot(filament.xyz_normal[k],filament.xyz_faces[k,2,:]) / bottom
		# vector from corner point (vertice 0) to point of intersection
		vect_corner_to_intersection	= radii_intersect[k]*sky.r_unit_vectors[idx_pixel] - filament.xyz_faces[k,0]
		projection10	= np.dot(vect_corner_to_intersection,filament.xyz_edges_vectors_unit[k,0])
		projection30	= np.dot(vect_corner_to_intersection,filament.xyz_edges_vectors_unit[k,1])
		if 0<=projection10<=np.linalg.norm(filament.xyz_edges_vectors[k,0,:]) and 0<=projection30<=np.linalg.norm(filament.xyz_edges_vectors[k,1,:]):
			# face k does intersect the ray r*hat(r)
			id_faces.append(k)
	#print(id_faces)
	return radii_intersect[id_faces]

def density(r,filament,sky,idx_pixel):
	# we need to transform r*hat(r) to the XYZ coordinate system
	XYZ_coord	= np.matmul(filament.inv_rot_matrix,r*sky.r_unit_vectors[idx_pixel] - np.array(filament.center))
	radius		= np.sqrt((XYZ_coord[0]/filament.sizes[0])**2 + (XYZ_coord[1]/filament.sizes[1])**2 + (XYZ_coord[2]/filament.sizes[2])**2)
	profile		= np.exp(-radius)
	return profile

def Bxyz(r_array,filament,sky,idx_pixel,magfield):
	# Project the local magnetic field into the local triad
	Nr					= len(r_array)
	result 				= np.zeros((Nr,4))
	local_magfield		= np.array([magfield.interp_fn((r*sky.r_unit_vectors[idx_pixel,0],r*sky.r_unit_vectors[idx_pixel,1],r*sky.r_unit_vectors[idx_pixel,2])) for r in r_array])
	# Bx,By,Bz,norm2
	result[:,0]			= np.dot(local_magfield,sky.local_triad[idx_pixel,0,:])
	result[:,1]			= np.dot(local_magfield,sky.local_triad[idx_pixel,1,:])
	result[:,2]			= np.dot(local_magfield,sky.local_triad[idx_pixel,2,:])
	result[:,3]			= result[:,0]**2 + result[:,1]**2 + result[:,2]**2
	return result

def integrator(filament,sky,idx_pixel,magfield,r1,r2,n=5):
	# an implementation of fixed_quad
	if r1 > r2:
		a=r2
		b=r1
	elif r2>r1:
		a=r1
		b=r2
	# This will integrate along the filament
	[x,w] 			= p_roots(n)
	x 				= np.real(x)
	y 				= (b-a)*(x+1)/2.0 + a
	density_value	= np.array([density(yy,filament,sky,idx_pixel) for yy in y])
	result			= Bxyz(y,filament,sky,idx_pixel,magfield)
	func_t			= density_value
	func_q			= density_value*(result[:,1]**2 - result[:,0]**2)/result[:,3]
	func_u			= density_value*(-2)*result[:,0]*result[:,1]/result[:,3]
	int_t			= (b-a)/2.0*sum(w*func_t,0)
	int_q			= (b-a)/2.0*sum(w*func_q,0)
	int_u			= (b-a)/2.0*sum(w*func_u,0)	
	return int_t,int_q,int_u

def paint_filament(n,sky,population,magfield):
	sizes 			= population.sizes[n]
	angles 			= population.angles[n]
	center 			= population.centers[n]
	filament		= Filament(center,sizes,angles)
#	print(filament.xyz_normal)
	pix_filament	= filament.do_query_polygon(sky)
	T_map			= np.zeros(12*sky.nside**2)
	Q_map			= np.zeros(12*sky.nside**2)
	U_map			= np.zeros(12*sky.nside**2)
	# Calculation of radii of intersection for each of the pixels in pix_filament
	N_pix_filament	= len(pix_filament)
	print('Filament %i has %i pixels'%(n,N_pix_filament))
	for n in range(N_pix_filament):
		idx_pixel			= pix_filament[n]
		# r1 and r2 are the radial distances from the origin to the 2 faces that intersect the cuboid
		r_distances			= distances(filament,sky,idx_pixel)
		# check that they are 2
		if len(r_distances) != 2:
			print('Error: something went wrong in the intersection of rays and faces')
			exit()
		r1	= r_distances[0]
		r2	= r_distances[1]		
		int_t,int_q,int_u		= integrator(filament,sky,idx_pixel,magfield,r1,r2,5)
		T_map[idx_pixel]		+= int_t
		Q_map[idx_pixel]		+= int_q
		U_map[idx_pixel]		+= int_u
	return T_map,Q_map,U_map

def paint_filament_aux(args):
	return paint_filament(*args)
