#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <FilamentPaint.h>

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

double** FilamentPaint_xyzVertices(double rot_matrix[3][3], double sizes_arr[3], double centers_arr[3]){
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

