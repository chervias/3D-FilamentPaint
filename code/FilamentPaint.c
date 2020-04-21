#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <query_polygon_wrapper.h>

double** FilamentPaint_RotationMatrix(double angles_arr[2]){
	static double rot_matrix[3][3];
	double ca = cos(angles_arr[0]); 
	double cb = cos(angles_arr[1]);
	double sa = sin(angles_arr[0]); 
	double sb = sin(angles_arr[1]);
	rot_matrix[0][0] = ca*cb;
	rot_matrix[1][0] = sa*cb;
	rot_matrix[2][0] = -sb;
	rot_matrix[0][1] = -sa;
	rot_matrix[1][1] = ca;
	rot_matrix[2][1] = 0;
	rot_matrix[0][2] = ca*sb;
	rot_matrix[1][2] = sa*sb;
	rot_matrix[2][2] = cb;
	return rot_matrix;
}
double** FilamentPaint_InvertRotMat(double rot_matrix[3][3]){
	// Following wikipedia https://en.wikipedia.org/wiki/Invertible_matrix#Methods_of_matrix_inversion
	static double inv_rot_matrix[3][3];
	double determinant=0.0;
	int i,j;
	for (i=0;i<3;i++){
		determinant = determinant + (rot_matrix[0][i]*(rot_matrix[1][(i+1)%3]*rot_matrix[2][(i+2)%3] - rot_matrix[1][(i+2)%3]*rot_matrix[2][(i+1)%3]));
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			inv_rot_matrix[j][i] = ((rot_matrix[(i+1)%3][(j+1)%3] * rot_matrix[(i+2)%3][(j+2)%3]) - (rot_matrix[(i+1)%3][(j+2)%3]*rot_matrix[(i+2)%3][(j+1)%3]))/ determinant;
		}
	}
	return inv_rot_matrix;
}

double** FilamentPaint_MatMul(double matrix_a[3][3],double matrix_b[3][3]){
	static double mat[3][3];
	int i,j,k;
	for (i=0;i<3;i++){
		for (j=0;j<3;j++){
			double sum = 0.0 ;
			for (k=0;k<3;k++){
				sum = sum + matrix_a[i][k]*matrix_b[k][j];
			}
			printf("%d \n",sum) ;
			mat[i][j] = sum ;
		}
	}
	return mat ;
}

double** FilamentPaint_xyzVertices(double rot_matrix[3][3], double sizes_arr[3], double centers_arr[3], double Size, int* isInside){
	/* This calculates the vertices of the cuboid in the xyz fixed coord */
	static double xyz_vertices[8][3],XYZ_vertices[8][3];
	int i,j,k;
	XYZ_vertices[0][0] = XYZ_vertices[1][0] = XYZ_vertices[2][0]=XYZ_vertices[3][0] = +5*sizes_arr[0];
	/* X element of vertices 5,6,7,8 */
	XYZ_vertices[4][0] = XYZ_vertices[5][0] = XYZ_vertices[6][0]=XYZ_vertices[7][0] = -5*sizes_arr[0];
	/* Y element of vertices 1,4,5,8 */
	XYZ_vertices[0][1] = XYZ_vertices[3][1] = XYZ_vertices[4][1]=XYZ_vertices[7][1] = -5*sizes_arr[1];
	/* Y elemeny of vertices 2,3,6,7 */
	XYZ_vertices[1][1] = XYZ_vertices[2][1] = XYZ_vertices[5][1]=XYZ_vertices[6][1] = +5*sizes_arr[1];
	/* Z element of vertices 1,2,5,6 */
	XYZ_vertices[0][2] = XYZ_vertices[1][2] = XYZ_vertices[4][2]=XYZ_vertices[5][2] = +5*sizes_arr[2];
	/* Z element of vertices 3,4,7,8 */
	XYZ_vertices[2][2] = XYZ_vertices[3][2] = XYZ_vertices[6][2]=XYZ_vertices[7][2] = -5*sizes_arr[2];
	
	/* multiply rot matrix by XYZ vertex and move by center_arr */
	for (i=0;i<8;i++){
		// j is rows
		for (j=0;j<3;j++){
			double sum = 0.0;
			for (k=0;k<3;k++){
				sum = sum + rot_matrix[j][k]*XYZ_vertices[i][k] ;
			}
			xyz_vertices[i][j] = sum + centers_arr[j];
		}
	}
	
	// Check if the vertices are outside the box
	// It will be 1. It will be changed to 0 if the condition is True
	for (i=0;i<8;i++){
		for (j=0;j<3;j++){
			if (xyz_vertices[i][j] < -0.5*Size || +0.5*Size < xyz_vertices[i][j]){
				*isInside = 0 ;
			}
		}
	}
	return xyz_vertices;
}

double** FilamentPaint_xyzNormalToFaces(double rot_matrix[3][3]){
	static double xyz_normal_to_faces[6][3],XYZ_normals[6][3];
	int i,j,k;
	XYZ_normals[0][0] = 1.0 ; XYZ_normals[0][1] = 0.0 ; XYZ_normals[0][2] = 0.0 ;
	XYZ_normals[1][0] = 1.0 ; XYZ_normals[1][1] = 0.0 ; XYZ_normals[1][2] = 0.0 ;
	XYZ_normals[2][0] = 0.0 ; XYZ_normals[2][1] = 1.0 ; XYZ_normals[2][2] = 0.0 ;
	XYZ_normals[3][0] = 0.0 ; XYZ_normals[3][1] = 1.0 ; XYZ_normals[3][2] = 0.0 ;
	XYZ_normals[4][0] = 0.0 ; XYZ_normals[4][1] = 0.0 ; XYZ_normals[4][2] = 1.0 ;
	XYZ_normals[5][0] = 0.0 ; XYZ_normals[5][1] = 0.0 ; XYZ_normals[5][2] = 1.0 ;
	
	/* multiply rot matrix by XYZ normals */
	for (i=0;i<6;i++){
		// j is rows
		for (j=0;j<3;j++){
			double sum = 0.0;
			for (k=0;k<3;k++){
				sum = sum + rot_matrix[j][k]*XYZ_normals[i][k] ;
			}
			xyz_normal_to_faces[i][j] = sum;
		}
	}
	return xyz_normal_to_faces;
}

double*** FilamentPaint_xyzFaces(double xyz_vertices[8][3]){
	/* The faces are 6, and each one has 4 vertices, each vertex has x,y,z*/
	static double xyz_faces[6][4][3];
	int i;
	for (i=0;i<3;i++){
		/* face 0 has vertices 0,1,2,3 */
		xyz_faces[0][0][i]	= xyz_vertices[0][i];
		xyz_faces[0][1][i]	= xyz_vertices[1][i];
		xyz_faces[0][2][i]	= xyz_vertices[2][i];
		xyz_faces[0][3][i]	= xyz_vertices[3][i];
		/* face 1 has vertices 4,5,6,7 */
		xyz_faces[1][0][i]	= xyz_vertices[4][i];
		xyz_faces[1][1][i]	= xyz_vertices[5][i];
		xyz_faces[1][2][i]	= xyz_vertices[6][i];
		xyz_faces[1][3][i]	= xyz_vertices[7][i];
		/* face 2 has vertices 0,3,7,4 */
		xyz_faces[2][0][i]	= xyz_vertices[0][i];
		xyz_faces[2][1][i]	= xyz_vertices[3][i];
		xyz_faces[2][2][i]	= xyz_vertices[7][i];
		xyz_faces[2][3][i]	= xyz_vertices[4][i];
		/* face 3 has vertices 1,2,6,5 */
		xyz_faces[3][0][i]	= xyz_vertices[1][i];
		xyz_faces[3][1][i]	= xyz_vertices[2][i];
		xyz_faces[3][2][i]	= xyz_vertices[6][i];
		xyz_faces[3][3][i]	= xyz_vertices[5][i];
		/* face 4 has vertices 0,1,5,4 */
		xyz_faces[4][0][i]	= xyz_vertices[0][i];
		xyz_faces[4][1][i]	= xyz_vertices[1][i];
		xyz_faces[4][2][i]	= xyz_vertices[5][i];
		xyz_faces[4][3][i]	= xyz_vertices[4][i];
		/* face 5 has vertices 2,3,7,6 */
		xyz_faces[5][0][i]	= xyz_vertices[2][i];
		xyz_faces[5][1][i]	= xyz_vertices[3][i];
		xyz_faces[5][2][i]	= xyz_vertices[7][i];
		xyz_faces[5][3][i]	= xyz_vertices[6][i];
	}
	return xyz_faces;
}

double*** FilamentPaint_xyzEdgeVectors(double xyz_faces[6][4][3]){
	static double xyz_edges[6][2][3];
	int i,k;
	for (i=0;i<6;i++){
		for (k=0;k<3;k++){
			xyz_edges[i][0][k]	= xyz_faces[i][1][k] - xyz_faces[i][0][k];
			xyz_edges[i][1][k]	= xyz_faces[i][3][k] - xyz_faces[i][0][k];
		}
	}
	return xyz_edges;
}
double*** FilamentPaint_xyzEdgeVectorsUnit(double xyz_edges[6][2][3]){
	static double xyz_edges_unit[6][2][3];
	int i,j ;
	for (i=0;i<6;i++){
		double norm0=0.0,norm1=0.0;
		for (j=0;j<3;j++){
			norm0 = norm0 + pow(xyz_edges[i][0][j],2);
			norm1 = norm1 + pow(xyz_edges[i][1][j],2);
		}
		for (j=0;j<3;j++){
			xyz_edges_unit[i][0][j] = xyz_edges[i][0][j]/sqrt(norm0);
			xyz_edges_unit[i][1][j] = xyz_edges[i][1][j]/sqrt(norm1);
		}
	}
	return xyz_edges_unit;
}
void FilamentPaint_DoQueryPolygon(long nside_long,double xyz_faces[6][4][3], long* ipix, long* nipix){
	int i,j,k;
	int npix = 12*nside_long*nside_long;
    // This buffer is to keep the result from the 6 runs of query_polygon
    long* ipix_buffer = calloc(npix,sizeof(long));
    long nipix_buffer ;
    // nipix_current stores the current length of ipix instantaniously
    // nipix_previous stores the length of ipix up to the previous face
    int nipix_current=0,nipix_previous=0;
    for (i=0;i<6;i++){
        // transform the vector to an angle
        static double theta[4],phi[4];
        for (j=0;j<4;j++){
            double radius=0.0;
            for (k=0;k<3;k++){
                radius = radius + pow(xyz_faces[i][j][k],2);
            }
            phi[j]	= atan2(xyz_faces[i][j][1],xyz_faces[i][j][0]);
            theta[j]		= acos(xyz_faces[i][j][2]/sqrt(radius));
        }
        query_polygon_wrapper(theta,phi,nside_long,ipix_buffer,&nipix_buffer);
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

double* FilamentPaint_CalculateDistances(double xyz_normal_to_faces[6][3], double xyz_faces[6][4][3], double xyz_edges[6][2][3], double xyz_edges_unit[6][2][3], double rUnitVector_ipix[3]){
	// This function receives the Filament matrices and the r_unit_vector of the corresponding pixel and returns the 2 distances of intersection
	static double radii_intersect[6];
	static int id_faces[2];
	int i,j,c=0;
	for (i=0 ; i<6 ; i++){
		// iterate over the 6 faces
		// "bottom" is the dot product between the normal of face i and r_unit_vector, "top" is the dot product between the normal and some point in the plane
		double bottom=0.0,top=0.0;
		for (j=0 ; j<3 ; j++){bottom = bottom + xyz_normal_to_faces[i][j]*rUnitVector_ipix[j];}
		if (bottom==0.0){
			// if bottom is 0.0, the radius will never intersect the plane
			radii_intersect[i] = -1.0;
		}
		else{
			for (j=0;j<3;j++){top = top + xyz_normal_to_faces[i][j]*xyz_faces[i][2][j] ;}
			radii_intersect[i] = top / bottom ;
			static double VectCornerToIntersection[3];
			for (j=0;j<3;j++){
				VectCornerToIntersection[j] = radii_intersect[i]*rUnitVector_ipix[j] - xyz_faces[i][0][j] ;
			}
			double proj10 = 0.0, proj30=0.0 , norm0=0.0 , norm1=0.0;
			for (j=0;j<3;j++){
				proj10 = proj10 + VectCornerToIntersection[j] * xyz_edges_unit[i][0][j] ;
				proj30 = proj30 + VectCornerToIntersection[j] * xyz_edges_unit[i][1][j] ;
				norm0  = norm0 + pow(xyz_edges[i][0][j],2);
				norm1  = norm1 + pow(xyz_edges[i][1][j],2);
			}
			if ((0.0<=proj10) && (proj10<=sqrt(norm0)) && (0.0<=proj30) && (proj30<=sqrt(norm1))){
				// face i does intersect the ray r*hat(r)
				id_faces[c] = i ;
				c++;
			}
		}
	}
	static double distances[2];
	distances[0] = radii_intersect[id_faces[0]];
	distances[1] = radii_intersect[id_faces[1]];
	return distances ;
}

double FilamentPaint_Density(double r, double rot_matrix[3][3], double rUnitVector_ipix[3], double centers_arr[3], double sizes_arr[3]){
	// This function defines the density profile within the cuboid in the XYZ coordinates
	// rot_matrix must be the inverse rotation matrix
	int i,j;
	static double XYZ_coord[3];
	double radius, profile;
	for (i=0;i<3;i++){
		double sum=0.0 ;
		for (j=0;j<3;j++){
			sum = sum + rot_matrix[i][j]*(r*rUnitVector_ipix[j]-centers_arr[j]) ;
		}
		XYZ_coord[i] = sum ;
	}
	radius		= pow(XYZ_coord[0]/sizes_arr[0],2)+pow(XYZ_coord[1]/sizes_arr[1],2)+pow(XYZ_coord[2]/sizes_arr[2],2);
	profile 	= exp(-sqrt(radius)) ;
	return profile;
}

double* FilamentPaint_TrilinearInterpolation(PyObject* Bcube_obj, double size_box, int nbox, double vector[3]){
	// Bcube is the cube with the values, this cube has dimensions nbox*nbox*nbox pixels, size_box is the physical size of each side
	// vector is the vector for which we want to know the interpolated value
	// This follows https://en.wikipedia.org/wiki/Trilinear_interpolation
	int i;
	// First, we need to get the indices of the cube where vector lives
	int idx_x1 	= ceil(vector[0]*(nbox-1)/size_box + 0.5*(nbox-1.0));
	int idx_x0 	= floor(vector[0]*(nbox-1)/size_box + 0.5*(nbox-1.0));
	int idx_y1 	= ceil(vector[1]*(nbox-1)/size_box + 0.5*(nbox-1.0));
	int idx_y0 	= floor(vector[1]*(nbox-1)/size_box + 0.5*(nbox-1.0));
	int idx_z1 	= ceil(vector[2]*(nbox-1)/size_box + 0.5*(nbox-1.0));
	int idx_z0 	= floor(vector[2]*(nbox-1)/size_box + 0.5*(nbox-1.0));

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
	static double c[3];
	for (i=0;i<3;i++){
		c[i]		= c0[i]*(1.0-zd) + c1[i]*zd ;
	}
	return c;
}
double* FilamentPaint_Blocal(PyObject* Bcube_obj,double size_box, int nbox,double center_arr[3]){
	double* localB	= FilamentPaint_TrilinearInterpolation(Bcube_obj,size_box,nbox,center_arr);
	static double result[3];
	result[0] = localB[0];
	result[1] = localB[1];
	result[2] = localB[2];
	return result;
}

double* FilamentPaint_Bxyz(double r, PyObject* Bcube_obj,double size_box, int nbox,double rUnitVector_ipix[3],double local_triad[3][3]){
	// Get the local magnetic field in r*hat(r)
	int i;
	double vec[3] ;
	for (i=0 ; i<3 ; i++)
		vec[i]	= r*rUnitVector_ipix[i] ;
	double* localB	= FilamentPaint_TrilinearInterpolation(Bcube_obj,size_box,nbox,vec);
	// The result is a 1D array with size 4: Bx, By, Bz, norm2
	static double result[4];
	double Bx=0.0,By=0.0,Bz=0.0;
	for (i=0;i<3;i++){
		Bx	= Bx + localB[i]*local_triad[0][i] ;
		By	= By + localB[i]*local_triad[1][i] ;
		Bz	= Bz + localB[i]*local_triad[2][i] ;
	}
	result[0]	= Bx;
	result[1]	= By; 
	result[2]	= Bz;
	result[3]	= Bx*Bx + By*By + Bz*Bz;
	return result;
}
double* FilamentPaint_Integrator(double r1, double r2, double rot_matrix[3][3], double rUnitVector_ipix[3],double local_triad[3][3], double centers_arr[3], double sizes_arr[3], PyObject* Bcube_obj,double size_box, int nbox){
	// Integrator
	static double integ[3] ;
	double a,b ;
	int i;
	if (r1 > r2){a=r2;b=r1;}
	else{a=r1;b=r2;}
	// This is an implementation of python's fixed_quad
	// For n=5
	double x[5]		= {-9.061798459386639637003213465505E-01,-5.384693101056831077144693153969E-01,0.0E+00,5.384693101056831077144693153969E-01,9.061798459386639637003213465505E-01} ;
	double w[5]		= {2.369268850561889738770560143166E-01,4.786286704993665264140645376756E-01,5.688888888888891104400613585312E-01,4.786286704993665264140645376756E-01,2.369268850561889738770560143166E-01};
	double sumT=0.0, sumQ=0.0, sumU=0.0 ;
	for (i=0;i<5;i++){
		double y		= (b-a)*(x[i]+1)/2.0 + a ;
		// this now includes Larson's law
		double density_0 = pow(sizes_arr[2],-1.1) ;
		double density	= FilamentPaint_Density(y,rot_matrix,rUnitVector_ipix,centers_arr,sizes_arr) ;
		double* result	= FilamentPaint_Bxyz(y,Bcube_obj,size_box,nbox,rUnitVector_ipix,local_triad) ;
		sumT 			= sumT + w[i]*density_0*density ;
		sumQ			= sumQ + w[i]*density_0*density*(pow(result[1],2) - pow(result[0],2))/result[3] ;
		sumU			= sumU + w[i]*density_0*density*(-2.0)*result[1]*result[0]/result[3];
	}
	integ[0] 		= (b-a)/2.0*sumT ;
	integ[1]		= (b-a)/2.0*sumQ ;
	integ[2]		= (b-a)/2.0*sumU ;
	return integ;
}
