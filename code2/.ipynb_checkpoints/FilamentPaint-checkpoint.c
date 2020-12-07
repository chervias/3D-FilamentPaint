#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <query_polygon_wrapper.h>

void FilamentPaint_RotationMatrix(double* angles_arr, double* rot_matrix){
	double ca = cos(angles_arr[0]); 
	double cb = cos(angles_arr[1]);
	double sa = sin(angles_arr[0]); 
	double sb = sin(angles_arr[1]);
	rot_matrix[0*3+0] = ca*cb;
	rot_matrix[1*3+0] = sa*cb;
	rot_matrix[2*3+0] = -sb;
	rot_matrix[0*3+1] = -sa;
	rot_matrix[1*3+1] = ca;
	rot_matrix[2*3+1] = 0;
	rot_matrix[0*3+2] = ca*sb;
	rot_matrix[1*3+2] = sa*sb;
	rot_matrix[2*3+2] = cb;
}
void FilamentPaint_InvertRotMat(double* rot_matrix, double* inv_rot_matrix){
	// Following wikipedia https://en.wikipedia.org/wiki/Invertible_matrix#Methods_of_matrix_inversion
	double determinant=0.0;
	int i,j;
	for (i=0;i<3;i++){
		determinant = determinant + (rot_matrix[0*3+i]*(rot_matrix[1*3+(i+1)%3]*rot_matrix[2*3+(i+2)%3] - rot_matrix[1*3+(i+2)%3]*rot_matrix[2*3+(i+1)%3]));
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			inv_rot_matrix[j*3+i] = ((rot_matrix[((i+1)%3)*3+(j+1)%3] * rot_matrix[((i+2)%3)*3+(j+2)%3]) - (rot_matrix[((i+1)%3)*3+(j+2)%3]*rot_matrix[((i+2)%3)*3+(j+1)%3]))/ determinant;
		}
	}
}

void FilamentPaint_xyzVertices(double* rot_matrix, double* sizes_arr, double* centers_arr, double Size, unsigned int* isInside,double* xyz_vertices){
	/* This calculates the vertices of the cuboid in the xyz fixed coord */
	unsigned int i,j,k;
	double XYZ_vertices[24];
	XYZ_vertices[0*3+0] = XYZ_vertices[1*3+0] = XYZ_vertices[2*3+0]=XYZ_vertices[3*3+0] = +5*sizes_arr[0];
	/* X element of vertices 5,6,7,8 */
	XYZ_vertices[4*3+0] = XYZ_vertices[5*3+0] = XYZ_vertices[6*3+0]=XYZ_vertices[7*3+0] = -5*sizes_arr[0];
	/* Y element of vertices 1,4,5,8 */
	XYZ_vertices[0*3+1] = XYZ_vertices[3*3+1] = XYZ_vertices[4*3+1]=XYZ_vertices[7*3+1] = -5*sizes_arr[1];
	/* Y elemeny of vertices 2,3,6,7 */
	XYZ_vertices[1*3+1] = XYZ_vertices[2*3+1] = XYZ_vertices[5*3+1]=XYZ_vertices[6*3+1] = +5*sizes_arr[1];
	/* Z element of vertices 1,2,5,6 */
	XYZ_vertices[0*3+2] = XYZ_vertices[1*3+2] = XYZ_vertices[4*3+2]=XYZ_vertices[5*3+2] = +5*sizes_arr[2];
	/* Z element of vertices 3,4,7,8 */
	XYZ_vertices[2*3+2] = XYZ_vertices[3*3+2] = XYZ_vertices[6*3+2]=XYZ_vertices[7*3+2] = -5*sizes_arr[2];
	
	/* multiply rot matrix by XYZ vertex and move by center_arr */
	for (i=0;i<8;i++){
		// j is rows
		for (j=0;j<3;j++){
			double sum = 0.0;
			for (k=0;k<3;k++){
				sum = sum + rot_matrix[j*3+k]*XYZ_vertices[i*3+k] ;
			}
			xyz_vertices[i*3+j] = sum + centers_arr[j];
		}
	}
	
	// Check if the vertices are outside the box
	// It will be 1. It will be changed to 0 if the condition is True
	for (i=0;i<8;i++){
		for (j=0;j<3;j++){
			if (xyz_vertices[i*3+j] < -0.5*Size || +0.5*Size < xyz_vertices[i*3+j]){
				*isInside = 0 ;
			}
		}
	}
}

void FilamentPaint_xyzNormalToFaces(double* rot_matrix, double* xyz_normal_to_faces){
	//xyz_normal_to_faces has shape [6][3]
	double XYZ_normals[6*3];
	unsigned int i,j,k;
	XYZ_normals[0*3+0] = 1.0 ; XYZ_normals[0*3+1] = 0.0 ; XYZ_normals[0*3+2] = 0.0 ;
	XYZ_normals[1*3+0] = 1.0 ; XYZ_normals[1*3+1] = 0.0 ; XYZ_normals[1*3+2] = 0.0 ;
	XYZ_normals[2*3+0] = 0.0 ; XYZ_normals[2*3+1] = 1.0 ; XYZ_normals[2*3+2] = 0.0 ;
	XYZ_normals[3*3+0] = 0.0 ; XYZ_normals[3*3+1] = 1.0 ; XYZ_normals[3*3+2] = 0.0 ;
	XYZ_normals[4*3+0] = 0.0 ; XYZ_normals[4*3+1] = 0.0 ; XYZ_normals[4*3+2] = 1.0 ;
	XYZ_normals[5*3+0] = 0.0 ; XYZ_normals[5*3+1] = 0.0 ; XYZ_normals[5*3+2] = 1.0 ;
	
	/* multiply rot matrix by XYZ normals */
	for (i=0;i<6;i++){
		// j is rows
		for (j=0;j<3;j++){
			double sum = 0.0;
			for (k=0;k<3;k++){
				sum = sum + rot_matrix[j*3+k]*XYZ_normals[i*3+k] ;
			}
			xyz_normal_to_faces[i*3+j] = sum;
		}
	}
}

void FilamentPaint_xyzFaces(double* xyz_vertices, double* xyz_faces){
	/* The faces are 6, and each one has 4 vertices, each vertex has x,y,z*/
	//static double xyz_faces[6][4][3];
	// xyz_vertices is [8][3]
	unsigned int i;
	for (i=0;i<3;i++){
		/* face 0 has vertices 0,1,2,3 */
		xyz_faces[0*4*3+0*3+i]	= xyz_vertices[0*3+i];
		xyz_faces[0*4*3+1*3+i]	= xyz_vertices[1*3+i];
		xyz_faces[0*4*3+2*3+i]	= xyz_vertices[2*3+i];
		xyz_faces[0*4*3+3*3+i]	= xyz_vertices[3*3+i];
		/* face 1 has vertices 4,5,6,7 */
		xyz_faces[1*4*3+0*3+i]	= xyz_vertices[4*3+i];
		xyz_faces[1*4*3+1*3+i]	= xyz_vertices[5*3+i];
		xyz_faces[1*4*3+2*3+i]	= xyz_vertices[6*3+i];
		xyz_faces[1*4*3+3*3+i]	= xyz_vertices[7*3+i];
		/* face 2 has vertices 0,3,7,4 */
		xyz_faces[2*4*3+0*3+i]	= xyz_vertices[0*3+i];
		xyz_faces[2*4*3+1*3+i]	= xyz_vertices[3*3+i];
		xyz_faces[2*4*3+2*3+i]	= xyz_vertices[7*3+i];
		xyz_faces[2*4*3+3*3+i]	= xyz_vertices[4*3+i];
		/* face 3 has vertices 1,2,6,5 */
		xyz_faces[3*4*3+0*3+i]	= xyz_vertices[1*3+i];
		xyz_faces[3*4*3+1*3+i]	= xyz_vertices[2*3+i];
		xyz_faces[3*4*3+2*3+i]	= xyz_vertices[6*3+i];
		xyz_faces[3*4*3+3*3+i]	= xyz_vertices[5*3+i];
		/* face 4 has vertices 0,1,5,4 */
		xyz_faces[4*4*3+0*3+i]	= xyz_vertices[0*3+i];
		xyz_faces[4*4*3+1*3+i]	= xyz_vertices[1*3+i];
		xyz_faces[4*4*3+2*3+i]	= xyz_vertices[5*3+i];
		xyz_faces[4*4*3+3*3+i]	= xyz_vertices[4*3+i];
		/* face 5 has vertices 2,3,7,6 */
		xyz_faces[5*4*3+0*3+i]	= xyz_vertices[2*3+i];
		xyz_faces[5*4*3+1*3+i]	= xyz_vertices[3*3+i];
		xyz_faces[5*4*3+2*3+i]	= xyz_vertices[7*3+i];
		xyz_faces[5*4*3+3*3+i]	= xyz_vertices[6*3+i];
	}
}

void FilamentPaint_xyzEdgeVectors(double* xyz_faces, double* xyz_edges){
	//xyz_faces has shape [6][4][3]
	//xyz_edges has shape [6][2][3]
	unsigned int i,k;
	for (i=0;i<6;i++){
		for (k=0;k<3;k++){
			xyz_edges[i*2*3+0*3+k]	= xyz_faces[i*4*3+1*3+k] - xyz_faces[i*4*3+0*3+k];
			xyz_edges[i*2*3+1*3+k]	= xyz_faces[i*4*3+3*3+k] - xyz_faces[i*4*3+0*3+k];
		}
	}
}

void FilamentPaint_xyzEdgeVectorsUnit(double* xyz_edges, double* xyz_edges_unit){
	// xyz_edges has shape [6][2][3]
	//xyz_edges_unit has shape [6][2][3];
	unsigned int i,j ;
	for (i=0;i<6;i++){
		double norm0=0.0,norm1=0.0;
		for (j=0;j<3;j++){
			norm0 = norm0 + pow(xyz_edges[i*2*3+0*3+j],2);
			norm1 = norm1 + pow(xyz_edges[i*2*3+1*3+j],2);
		}
		for (j=0;j<3;j++){
			xyz_edges_unit[i*2*3+0*3+j] = xyz_edges[i*2*3+0*3+j]/sqrt(norm0);
			xyz_edges_unit[i*2*3+1*3+j] = xyz_edges[i*2*3+1*3+j]/sqrt(norm1);
		}
	}
}
void FilamentPaint_DoQueryPolygon(unsigned int nside, double* xyz_faces, unsigned long* ipix, unsigned long* nipix){
	unsigned long i,j,k;
	unsigned int npix_max = 200000;
    // This buffer is to keep the result from the 6 runs of query_polygon
    unsigned long* ipix_buffer = calloc(npix_max,sizeof(long));
    unsigned long nipix_buffer ;
    // nipix_current stores the current length of ipix instantaniously
    // nipix_previous stores the length of ipix up to the previous face
    unsigned long nipix_current=0,nipix_previous=0;
	double test,a,b,c,s;
//	for (i=0;i<6;i++){
//		printf("Face %i \n",i);
//		for (j=0;j<4;j++){
//			printf("corner%i = (%.20E,%.20E,%.20E)\n",j,xyz_faces[i*4*3 + j*3 + 0],xyz_faces[i*4*3 + j*3 + 1],xyz_faces[i*4*3 + j*3 + 2]);
//		}
//	}
    for (i=0;i<6;i++){
        // transform the vector to an angle
        static double theta[4],phi[4];
        for (j=0;j<4;j++){
            double radius=0.0;
            for (k=0;k<3;k++){
                radius = radius + pow(xyz_faces[i*4*3 + j*3 + k],2);
            }
            phi[j]	= atan2(xyz_faces[i*4*3 + j*3 + 1],xyz_faces[i*4*3 + j*3 + 0]);
            theta[j]= acos(xyz_faces[i*4*3 + j*3 + 2]/sqrt(radius));
        }
		// check the area of the triangle by 3 corners of the face Heron's formula
		// this assumes x=theta and y=phi, flat sky
		a = sqrt( pow(theta[3] - theta[2],2) + pow(phi[3] - phi[2],2) );
		b = sqrt( pow(theta[2] - theta[1],2) + pow(phi[2] - phi[1],2) );
		c = sqrt( pow(theta[1] - theta[3],2) + pow(phi[1] - phi[3],2) );
		s = 0.5*(a+b+c);
		test = sqrt(s*(s-a)*(s-b)*(s-c)) ; 
		//printf("Area of the projected triangle %.6E \n",test);
		if(fabs(test) > 6e-12){query_polygon_wrapper(theta,phi,nside,ipix_buffer,&nipix_buffer);}
		else{
			//printf("skip this face\n");
			nipix_buffer = 0;
			continue;
		}
        if (i==0){
            // If we are in face 0, copy the contents of ipix_buffer into ipix
            for (j=0;j<nipix_buffer;j++){
                ipix[j] = ipix_buffer[j] ;
                nipix_current++;
            }
        }
        else{
            // If we are in face 1-5, we need to check that the elements are unique before copying into ipix
            for (j=0;j<nipix_buffer;j++){
                // we will check if ipix_buffer[j] is already in ipix
                int IsAlready = 0;
                for (k=0;k<nipix_previous;k++){
                    if (ipix_buffer[j] == ipix[k]){
                        // If this happens, the pixel index is already in ipix and we break, move on to the next
                        // element in ipix_buffer
                        IsAlready = 1;
                        break;
                    }
                }
                if (IsAlready==0){
                    ipix[nipix_current] = ipix_buffer[j];
                    nipix_current++;
                }
            }
        }
        // we need to update nipix_previous for the next face
        nipix_previous=nipix_current ;
    }
    // assign nipix
    *nipix =  (long) nipix_current ;
    // Free the memory
    free(ipix_buffer);
}

void FilamentPaint_DoLocalTriad(unsigned long* ipix_final, unsigned long nipix, unsigned int nside, double* localtriad){
	// localtriad is double* localtriad = calloc(nipix*3*3,sizeof(double)), which is created before
	localtriad_wrapper(ipix_final,nipix,nside,localtriad);
}

//void FilamentPaint_DoUpgradewithKernel(double fwhm, unsigned int nside_low, unsigned int nside_high, unsigned int npix_high, unsigned int nipix, unsigned int* ipix_final,double* integ, double* TQUmap){
//	upgrade_with_kernel(fwhm,nside_low,nside_high,npix_high,nipix,ipix_final,integ,TQUmap);
//}

void FilamentPaint_CalculateDistances(double* xyz_normal_to_faces, double* xyz_faces, double* xyz_edges, double* xyz_edges_unit, unsigned long idx_pix, double* local_triad, double* distances){
	// This function receives the Filament matrices and the local_triad and the pixel index and returns the 2 distances of intersection
	// shapes: xyz_normal_to_faces[6][3] , xyz_faces[6][4][3] , xyz_edges[6][2][3] , xyz_edges_unit[6][2][3]
	double radii_intersect[6];
	unsigned int id_faces[2];
	unsigned int i,j,c=0;
	for (i=0 ; i<6 ; i++){
		// iterate over the 6 faces
		// "bottom" is the dot product between the normal of face i and r_unit_vector, "top" is the dot product between the normal and some point in the plane
		double bottom=0.0,top=0.0;
		// here j indexes the x,y,z componentes
		for (j=0 ; j<3 ; j++) bottom = bottom + xyz_normal_to_faces[i*3+j]*local_triad[idx_pix*3*3 + 2*3 + j];
		if (bottom==0.0){
			// if bottom is 0.0, the radius will never intersect the plane
			radii_intersect[i] = -1.0;
		}
		else{
			for (j=0;j<3;j++) top = top + xyz_normal_to_faces[i*3+j]*xyz_faces[i*4*3 + 2*3 + j] ;
			radii_intersect[i] = top / bottom ;
			double VectCornerToIntersection[3];
			for (j=0;j<3;j++) VectCornerToIntersection[j] = radii_intersect[i]*local_triad[idx_pix*3*3 + 2*3 + j] - xyz_faces[i*4*3 + 0*3 + j] ;
			double proj10 = 0.0, proj30=0.0 , norm0=0.0 , norm1=0.0;
			for (j=0;j<3;j++){
				proj10 = proj10 + VectCornerToIntersection[j] * xyz_edges_unit[i*2*3 + 0*3 + j] ;
				proj30 = proj30 + VectCornerToIntersection[j] * xyz_edges_unit[i*2*3 + 1*3 + j] ;
				norm0  = norm0 + pow(xyz_edges[i*2*3 + 0*3 + j],2);
				norm1  = norm1 + pow(xyz_edges[i*2*3 + 1*3 + j],2);
			}
			if ((0.0<=proj10) && (proj10<=sqrt(norm0)) && (0.0<=proj30) && (proj30<=sqrt(norm1))){
				// face i does intersect the ray r*hat(r)
				id_faces[c] = i ;
				c++;
			}
		}
	}
	distances[0] = radii_intersect[id_faces[0]];
	distances[1] = radii_intersect[id_faces[1]];
}

double FilamentPaint_Density(double r, double* rot_matrix, unsigned long idx_pix, double* local_triad, double* centers_arr, double* sizes_arr){
	// This function defines the density profile within the cuboid in the XYZ coordinates
	// rot_matrix must be the inverse rotation matrix
	unsigned int i,j;
	double XYZ_coord[3];
	double radius, profile;
	for (i=0;i<3;i++){
		double sum=0.0 ;
		for (j=0;j<3;j++){
			sum = sum + rot_matrix[i*3 + j]*(r*local_triad[idx_pix*3*3 + 2*3 + j] - centers_arr[j]) ;
		}
		XYZ_coord[i] = sum ;
	}
	radius		= pow(XYZ_coord[0]/sizes_arr[0],2)+pow(XYZ_coord[1]/sizes_arr[1],2)+pow(XYZ_coord[2]/sizes_arr[2],2);
	//if (radius <= 1.0){
	//	profile 	= exp(-sqrt(radius)) ;
		profile 	= exp(-0.5*radius) ;
	//}
	//else{
	//	profile 	= 0.0;
	//}
	return profile;
}

void FilamentPaint_TrilinearInterpolation(PyObject* Bcube_obj, double size_box, unsigned int nbox, double* vector, double* c){
	// Bcube is the cube with the values, this cube has dimensions nbox*nbox*nbox pixels, size_box is the physical size of each side
	// vector is the vector for which we want to know the interpolated value
	// This follows https://en.wikipedia.org/wiki/Trilinear_interpolation
	unsigned int i;
	// First, we need to get the indices of the cube where vector lives
	unsigned int idx_x1 	= ceil(vector[0]*(nbox-1)/size_box + 0.5*(nbox-1.0));
	unsigned int idx_x0 	= floor(vector[0]*(nbox-1)/size_box + 0.5*(nbox-1.0));
	unsigned int idx_y1 	= ceil(vector[1]*(nbox-1)/size_box + 0.5*(nbox-1.0));
	unsigned int idx_y0 	= floor(vector[1]*(nbox-1)/size_box + 0.5*(nbox-1.0));
	unsigned int idx_z1 	= ceil(vector[2]*(nbox-1)/size_box + 0.5*(nbox-1.0));
	unsigned int idx_z0 	= floor(vector[2]*(nbox-1)/size_box + 0.5*(nbox-1.0));

	// map the indices to real coordinates
	double x0		= size_box*idx_x0/(nbox-1.) - 0.5*size_box ;
	double x1		= size_box*idx_x1/(nbox-1.) - 0.5*size_box ;
	double y0		= size_box*idx_y0/(nbox-1.) - 0.5*size_box ;
	double y1		= size_box*idx_y1/(nbox-1.0) - 0.5*size_box ;	
	double z0		= size_box*idx_z0/(nbox-1.0) - 0.5*size_box ;
	double z1		= size_box*idx_z1/(nbox-1.0) - 0.5*size_box ;
	
	// Calculate xd,yd,zd
	double xd		= (vector[0] - x0)/(x1 - x0) ;
	double yd		= (vector[1] - y0)/(y1 - y0) ;
	double zd		= (vector[2] - z0)/(z1 - z0) ;
	
	// interpolate along x
	double c00[3],c01[3],c10[3],c11[3] ;
	for (i=0;i<3;i++){
		c00[i]		= (*(double*)PyArray_GETPTR4(Bcube_obj,idx_x0,idx_y0,idx_z0,i))*(1.0 - xd) 
					+ (*(double*)PyArray_GETPTR4(Bcube_obj,idx_x1,idx_y0,idx_z0,i))*xd ;
		c01[i]		= (*(double*)PyArray_GETPTR4(Bcube_obj,idx_x0,idx_y0,idx_z1,i))*(1.0 - xd) 
					+ (*(double*)PyArray_GETPTR4(Bcube_obj,idx_x1,idx_y0,idx_z1,i))*xd ;
		c10[i]		= (*(double*)PyArray_GETPTR4(Bcube_obj,idx_x0,idx_y1,idx_z0,i))*(1.0 - xd) 
					+ (*(double*)PyArray_GETPTR4(Bcube_obj,idx_x1,idx_y1,idx_z0,i))*xd ;
		c11[i]		= (*(double*)PyArray_GETPTR4(Bcube_obj,idx_x0,idx_y1,idx_z1,i))*(1.0 - xd) 
					+ (*(double*)PyArray_GETPTR4(Bcube_obj,idx_x1,idx_y1,idx_z1,i))*xd ;
	}
	// interpolate along y
	double c0[3],c1[3];
	for (i=0;i<3;i++){
		c0[i]		= c00[i]*(1.0 - yd) + c10[i]*yd ;
		c1[i]		= c01[i]*(1.0 - yd) + c11[i]*yd ;
	}
	// interpolate along z
	for (i=0;i<3;i++){
		c[i]		= c0[i]*(1.0-zd) + c1[i]*zd ;
	}
}

void FilamentPaint_Bxyz(double r, PyObject* Bcube_obj,double size_box, unsigned int nbox, unsigned long idx_pix, double* local_triad, double* result){
	// Get the local magnetic field in r*hat(r)
	unsigned int i;
	double vec[3], localB[3] ;
	for (i=0 ; i<3 ; i++) vec[i] = r*local_triad[idx_pix*3*3 + 2*3 + i] ;
	FilamentPaint_TrilinearInterpolation(Bcube_obj,size_box,nbox,vec,localB);
	// The result is a 1D array with size 4: Bx, By, Bz, norm2
	double Bx=0.0,By=0.0,Bz=0.0;
	for (i=0;i<3;i++){
		Bx	= Bx + localB[i]*local_triad[idx_pix*3*3 + 0*3 + i] ;
		By	= By + localB[i]*local_triad[idx_pix*3*3 + 1*3 + i] ;
		Bz	= Bz + localB[i]*local_triad[idx_pix*3*3 + 2*3 + i] ;
	}
	result[0]	= Bx;
	result[1]	= By; 
	result[2]	= Bz;
	result[3]	= Bx*Bx + By*By + Bz*Bz;
}
void FilamentPaint_RiemannIntegrator(double r1, double r2, double* rot_matrix, unsigned long idx_pix, double* local_triad, double* centers_arr, double* sizes_arr, PyObject* Bcube_obj, double size_box, unsigned int nbox, double* integ, PyObject* fpol0, PyObject* thetaH){
	//Integrator
	double a,b ;
	unsigned int i;
	if (r1 > r2){a=r2;b=r1;}
	else{a=r1;b=r2;}
	double deltar = (b-a)/99.0;
	double sumT=0.0;
	double sumQ=0.0;
	double sumU=0.0;
	double fpol0_ = PyFloat_AsDouble(fpol0);
	double thetaH_ = PyFloat_AsDouble(thetaH);
	double fpol = fpol0_*pow(sin(thetaH_),2) ;
	double density_0 = pow(sizes_arr[2],-1.1);
	for (i=0;i<100;i++){
		double density	= FilamentPaint_Density(a+i*deltar,rot_matrix,idx_pix,local_triad,centers_arr,sizes_arr) ;
		double result[4];
		FilamentPaint_Bxyz(a+i*deltar,Bcube_obj,size_box,nbox,idx_pix,local_triad,result);
		sumT = sumT + density_0*density ;
		sumQ = sumQ + fpol*density_0*density*(pow(result[1],2) - pow(result[0],2))/result[3];
		sumU = sumU + fpol*density_0*density*(-2.0)*result[1]*result[0]/result[3];
	}
	integ[idx_pix*3 + 0] = sumT*deltar ;
	integ[idx_pix*3 + 1] = sumQ*deltar ;
	integ[idx_pix*3 + 2] = sumU*deltar ;
}
