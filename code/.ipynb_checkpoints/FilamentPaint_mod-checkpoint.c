#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <numpy/ndarrayobject.h>
#include <FilamentPaint.h>

static PyObject *Paint_Filament(PyObject *self, PyObject *args){
	/* Getting the elements */
	PyObject *Sizes_total=NULL,*Angles_total=NULL,*Centers_total=NULL,*Nside_value=NULL,*rUnitVectors=NULL ;
	PyObject *Npix_cube_obj=NULL,*Bcube_obj=NULL,*Size_obj=NULL,*LocalTriad_obj=NULL ;
	PyObject *n = NULL;
	PyObject *nside = NULL;
	PyObject *r_unit_vectors=NULL;
	PyObject *local_triad=NULL;
	PyObject *Population=NULL;
	PyObject *MagField=NULL;
	double sizes_arr[3], angles_arr[2], centers_arr[3] ;
	int i,j,k,isInside=1;
    
	if (!PyArg_ParseTuple(args, "OOOOOO",&n, &nside, &r_unit_vectors, &local_triad, &Population, &MagField))
		return NULL;
	
	/* Extract the attribute from Population object*/
	Sizes_total 		= PyObject_GetAttr(Population,PyUnicode_FromString("sizes"));	
	Angles_total 		= PyObject_GetAttr(Population,PyUnicode_FromString("angles"));
	Centers_total 		= PyObject_GetAttr(Population,PyUnicode_FromString("centers"));
	
	long nside_long     = PyLong_AsLong(nside);

	// Extract the attributes from MagField object
	Npix_cube_obj		= PyObject_GetAttr(MagField,PyUnicode_FromString("pixels")) ;
	int Npix_cube		= (int) PyLong_AsLong(Npix_cube_obj);
	Size_obj			= PyObject_GetAttr(MagField,PyUnicode_FromString("size")) ;
	double Size			= PyFloat_AsDouble(Size_obj);
	Bcube_obj			= PyObject_GetAttr(MagField,PyUnicode_FromString("Bcube")) ;
	/* We want the filament n*/
	long n_fil			= PyLong_AsLong(n);
	long npix			= 12*nside_long*nside_long;

	// Extract from Population object
	for (k=0;k<3;k++){
		sizes_arr[k]	= *(double*)PyArray_GETPTR2(Sizes_total, n_fil, k);
		centers_arr[k]	= *(double*)PyArray_GETPTR2(Centers_total, n_fil, k);
	}
	for (k=0;k<2;k++){
		angles_arr[k]	= *(double*)PyArray_GETPTR2(Angles_total, n_fil, k);
	}
	
	// initialize the healpix map
	double* TQUmap	= calloc(3*npix,sizeof(double));
	
	// This is for testing if the cuboid is outside the box
	//printf("is inside before = %i\n",isInside);	
	/* Calculate the rot matrix */
	double** rot_matrix 			= FilamentPaint_RotationMatrix(angles_arr);	
	double** inv_rot_matrix			= FilamentPaint_InvertRotMat(rot_matrix);
	/* Calculate the 8 vertices in the xyz coordinates */
	double** xyz_vertices			= FilamentPaint_xyzVertices(rot_matrix,sizes_arr,centers_arr,Size,&isInside);
	//printf("is inside after = %i\n",isInside);
	
	if (isInside==1){
		/* Calculate normal to faces rotated */
		double** xyz_normal_to_faces	= FilamentPaint_xyzNormalToFaces(rot_matrix);
		/* Calculate the faces matrix*/
		double*** xyz_faces				= FilamentPaint_xyzFaces(xyz_vertices);
		// Calculate the edges vectors
		double*** xyz_edges				= FilamentPaint_xyzEdgeVectors(xyz_faces);
		double*** xyz_edges_unit        = FilamentPaint_xyzEdgeVectorsUnit(xyz_edges);
		// Calculate the polygon 
		long* ipix = calloc(npix,sizeof(long)); 
		long nipix;
		FilamentPaint_DoQueryPolygon(nside_long,xyz_faces,ipix,&nipix);
		long* ipix_final = calloc(nipix,sizeof(long));
		for (j=0;j<nipix;j++)
			ipix_final[j] = ipix[j] ;
		free(ipix);
		
		// Cycle through each of the pixels in ipix
		printf("Filament %i has %i pixels \n",n_fil,nipix) ;
		for (i=0;i<nipix;i++){
			int index_pix = ipix_final[i];
			// Get the hat(r) vector for pixel index_pixel
			double rUnitVector_ipix[3];
			double LocalTriad_ipix[3][3];
			for (j=0;j<3;j++){
				rUnitVector_ipix[j] = *(double*)PyArray_GETPTR2(r_unit_vectors, index_pix, j);
				for (k=0;k<3;k++){
					LocalTriad_ipix[k][j]	= *(double*)PyArray_GETPTR3(local_triad,index_pix, k, j) ;
				}
			}
			double* rDistances 		= FilamentPaint_CalculateDistances(xyz_normal_to_faces,xyz_faces,xyz_edges,xyz_edges_unit,rUnitVector_ipix);
			double* integ			= FilamentPaint_Integrator(rDistances[0],rDistances[1],inv_rot_matrix,rUnitVector_ipix,LocalTriad_ipix,centers_arr,sizes_arr,Bcube_obj,Size,Npix_cube) ;

			for (j=0;j<3;j++){
				TQUmap[j*npix + index_pix]	= integ[j] ;
			}
		}
		free(ipix_final);
		
		npy_intp npy_shape[2] = {3,npix};
		PyObject *arr 		= PyArray_SimpleNewFromData(2,npy_shape, NPY_DOUBLE, TQUmap);
		PyArray_ENABLEFLAGS((PyArrayObject *)arr, NPY_OWNDATA);
		return(arr);
	}
	else{
		printf("Filament %i is outside the box, skipping\n",n_fil) ;
		npy_intp npy_shape[2] = {3,npix};
		PyObject *arr 		= PyArray_SimpleNewFromData(2,npy_shape, NPY_DOUBLE, TQUmap);
		PyArray_ENABLEFLAGS((PyArrayObject *)arr, NPY_OWNDATA);
		return(arr);
	}
}

static PyMethodDef FilamentPaintMethods[] = {
  {"Paint_Filament",  Paint_Filament, METH_VARARGS,NULL},
 {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef FilamentPaint_module = {
    PyModuleDef_HEAD_INIT,
    "FilamentPaint",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    FilamentPaintMethods
};

PyMODINIT_FUNC PyInit_FilamentPaint(void){
  PyObject *m;
  m = PyModule_Create(&FilamentPaint_module);
  import_array();  // This is important for using the numpy_array api, otherwise segfaults!
  return(m);
}
