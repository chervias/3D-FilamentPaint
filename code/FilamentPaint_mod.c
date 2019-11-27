#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <numpy/ndarrayobject.h>
#include <FilamentPaint.h>

static PyObject *Paint_Filament(PyObject *self, PyObject *args){
	/* Getting the elements */
	PyObject *sizes,*Sizes_total,*angles,*Angles_total,*centers,*Centers_total,*nside,*Nside_value;
	PyObject *n = NULL;
	PyObject *Sky=NULL;
	PyObject *Population=NULL;
	PyObject *MagField=NULL;
	double sizes_arr[3], angles_arr[2], centers_arr[3] ;
	int i,j,k;
	
	if (!PyArg_ParseTuple(args, "OOOO", &n, &Sky, &Population, &MagField))
		return NULL;
	/* Query for the attribute sizes of the Population object*/
	sizes				= PyUnicode_FromString("sizes");
	Sizes_total 		= PyObject_GetAttr(Population,sizes);	
	angles				= PyUnicode_FromString("angles");
	Angles_total 		= PyObject_GetAttr(Population,angles);
	centers				= PyUnicode_FromString("centers");
	Centers_total 		= PyObject_GetAttr(Population,centers);
	// Query the attributes on the Sky object
	nside				= PyUnicode_FromString("nside");
	Nside_value			= PyObject_GetAttr(Sky,nside);
	long nside_long		= PyLong_AsLong(Nside_value);
	
	
	/* We want the filament n*/
	long n_fil			= PyLong_AsLong(n);
	for (k=0;k<3;k++){
		sizes_arr[k]	= *(double*)PyArray_GETPTR2(Sizes_total, n_fil, k);
		centers_arr[k]	= *(double*)PyArray_GETPTR2(Centers_total, n_fil, k);
	}
	for (k=0;k<2;k++){
		angles_arr[k]	= *(double*)PyArray_GETPTR2(Angles_total, n_fil, k);
	}
	
	/* Calculate the rot matrix */
	double** rot_matrix 			= FilamentPaint_RotationMatrix(angles_arr);	
	/* Calculate the 8 vertices in the xyz coordinates */
	double** xyz_vertices			= FilamentPaint_xyzVertices(rot_matrix,sizes_arr,centers_arr);
	/* Calculate normal to faces rotated */
	double** xyz_normal_to_faces	= FilamentPaint_xyzNormalToFaces(rot_matrix);
	/* Calculate the faces matrix*/
	double*** xyz_faces				= FilamentPaint_xyzFaces(xyz_vertices);
	// Calculate the edges vectors
	double*** xyz_edges				= FilamentPaint_xyzEdgeVectors(xyz_faces);
	// Calculate the polygon 
	long *ipix						= calloc(12*nside_long*nside_long,sizeof(long));
	long n_ipix;
	ipix 							= FilamentPaint_DoQueryPolygon(nside_long,xyz_faces,n_ipix);
	
	printf("%d\n",n_ipix);
	
	int nd=3;
	npy_intp dims[3] = {6,2,3};
	PyObject *arr = PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE,xyz_edges);
	PyArray_ENABLEFLAGS((PyArrayObject *)arr, NPY_OWNDATA);	
	return(arr);
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

PyMODINIT_FUNC PyInit_FilamentPaint(void)
{
  PyObject *m;

  m = PyModule_Create(&FilamentPaint_module);
  
  import_array();  // This is important for using the numpy_array api, otherwise segfaults!

  return(m);
}
