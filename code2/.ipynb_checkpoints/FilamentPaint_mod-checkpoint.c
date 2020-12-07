#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <numpy/ndarrayobject.h>
#include <FilamentPaint.h>

static PyObject *Paint_Filament(PyObject *self, PyObject *args){
	/* Getting the elements */
	PyObject *Sizes_total=NULL,*Angles_total=NULL,*Centers_total=NULL ;
	PyObject *n = NULL;
	PyObject *nside = NULL;
	PyObject *Population=NULL;
	PyObject *Bcube=NULL;
	PyObject *size=NULL;
	PyObject *Npix_magfield=NULL;
	PyObject *resolution=NULL;
	PyObject *fpol0=NULL;
	PyObject *thetaH=NULL;
	double sizes_arr[3], angles_arr[2], centers_arr[3] ;
	unsigned long i;
	unsigned int j,k,isInside=1;
	// The inputs are n=filament number, nside, I don't need local triad anymore ? , Population object, Bcube, size, Npix_magfield
	if (!PyArg_ParseTuple(args, "OOOOOOOOO",&n, &nside, &Population, &Bcube, &size, &Npix_magfield,&resolution,&fpol0,&thetaH))
		return NULL;
	
	/* Extract the attribute from Population object*/
	Sizes_total 		= PyObject_GetAttr(Population,PyUnicode_FromString("sizes"));	
	Angles_total 		= PyObject_GetAttr(Population,PyUnicode_FromString("angles"));
	Centers_total 		= PyObject_GetAttr(Population,PyUnicode_FromString("centers"));
	
	unsigned int nside_fixed = (int) PyLong_AsLong(nside);
	unsigned long npix_fixed  = 12*nside_fixed*nside_fixed;
	unsigned int resol_int   = (int) PyLong_AsLong(resolution) ; 
	unsigned int Npix_magfield_int = (int) PyLong_AsLong(Npix_magfield);
	double Size_double    =  PyFloat_AsDouble(size);

	/* We want the filament n*/
	unsigned int n_fil = (int) PyLong_AsLong(n);

	// Extract from Population object
	for (k=0;k<3;k++){
		sizes_arr[k]	= *(double*)PyArray_GETPTR2(Sizes_total, n_fil, k);
		centers_arr[k]	= *(double*)PyArray_GETPTR2(Centers_total, n_fil, k);
	}
	for (k=0;k<2;k++){
		angles_arr[k]	= *(double*)PyArray_GETPTR2(Angles_total, n_fil, k);
	}
	
	// Calculate which resolution I need to sample in 50x50 pixels at least. 2^n_nside is the nside necesary for the sampling 
	unsigned int n_nside = round(log(0.1*resol_int*sqrt(M_PI/3.0)*sqrt(pow(centers_arr[0],2)+pow(centers_arr[1],2)+pow(centers_arr[2],2))/sizes_arr[0])/log(2.0)) ;
	unsigned int nside_variable = pow(2,n_nside);
	unsigned long npix_variable = 12*nside_variable*nside_variable ;

	// This is for testing if the cuboid is outside the box
	//printf("is inside before = %i\n",isInside);	
	/* Calculate the rot matrix */
	double rot_matrix[9],inv_rot_matrix[9],xyz_vertices[8*3], xyz_normal_to_faces[6*3], xyz_faces[6*4*3], xyz_edges[6*2*3], xyz_edges_unit[6*2*3];
	
	FilamentPaint_RotationMatrix(angles_arr,rot_matrix);
	FilamentPaint_InvertRotMat(rot_matrix,inv_rot_matrix);
	/* Calculate the 8 vertices in the xyz coordinates */
	FilamentPaint_xyzVertices(rot_matrix,sizes_arr,centers_arr,Size_double,&isInside,xyz_vertices);
	//printf("is inside after = %i\n",isInside);
	
	if (isInside==1){
		/* Calculate normal to faces rotated */
		FilamentPaint_xyzNormalToFaces(rot_matrix,xyz_normal_to_faces);
		/* Calculate the faces matrix*/
		FilamentPaint_xyzFaces(xyz_vertices,xyz_faces);
		// Calculate the edges vectors
		FilamentPaint_xyzEdgeVectors(xyz_faces,xyz_edges);
		FilamentPaint_xyzEdgeVectorsUnit(xyz_edges,xyz_edges_unit);
		//printf("Filament %i with nside %i \n",n_fil,nside_variable) ;
		// Calculate the polygon in the variable nside map
		unsigned long* ipix = calloc(200000,sizeof(long)); 
		unsigned long nipix;
		FilamentPaint_DoQueryPolygon(nside_variable,xyz_faces,ipix,&nipix);
		unsigned long* ipix_final = calloc(nipix,sizeof(long));
		for (j=0;j<nipix;j++)
			ipix_final[j] = ipix[j] ;
		free(ipix);
		// Cycle through each of the pixels in ipix
		printf("Filament %i has %lu pixels in a nside=%i pixelization \n",n_fil,nipix,nside_variable) ;
		double* localtriad = calloc(nipix*3*3,sizeof(double));
		FilamentPaint_DoLocalTriad(ipix_final,nipix,nside_variable,localtriad);
		// Now I won't use a full size map at nside_variable, because it would be a waste of memory. Instead, I will map the integ result into the nested nside_fixed map
		// integ will be shape [nipix][3]
		double* integ = calloc(nipix*3,sizeof(double));
		for (i=0;i<nipix;i++){
			double rDistances[2];
			FilamentPaint_CalculateDistances(xyz_normal_to_faces,xyz_faces,xyz_edges,xyz_edges_unit,i,localtriad,rDistances);
			FilamentPaint_RiemannIntegrator(rDistances[0],rDistances[1],inv_rot_matrix,i,localtriad,centers_arr,sizes_arr,Bcube,Size_double,Npix_magfield_int,integ,fpol0,thetaH);
		}
		// if we are upgrading, we have a special procedure
		// the map will be returned in the nside_variable resolution 
		if (nside_fixed > nside_variable){
			// This means the filament is bigger than the fixed resolution and we have to upgrade the map
			// we will have to cycle over all children pixels in the nside_fixed pixelization, assigning the same value to each
			// the first children pixel will be indexed by index_pix << 2*step, and the last children pixel will be indexed by (index_pix << 2*step)+pow(4,step)
			//unsigned int step = (int)(log(nside_fixed)/log(2.0) - log(nside_variable)/log(2.0)) ;
			//unsigned int index_children_pixel = index_pix << 2*step ;
			//for (k=index_children_pixel;k<(index_children_pixel+pow(4,step)+1);k++){
			//	for (j=0;j<3;j++) TQUmap[j*npix_fixed + k] = integ[j] ;
			//}
			// Now we will change how to do this
			// TQUmap will be modified 
			// fwhm in radians
			//double fwhm = sqrt(M_PI/3.0)/nside_variable ;
			//FilamentPaint_DoUpgradewithKernel(fwhm,nside_variable,nside_fixed,npix_fixed,nipix,ipix_final,integ,TQUmap) ;
			// in this case, TQUmap has a shape [3][npix_variable]
			double* TQUmap	= calloc(3*npix_variable,sizeof(double));
			for (i=0;i<nipix;i++){
				unsigned long index_pix = ipix_final[i];
				for (j=0;j<3;j++) TQUmap[j*npix_variable + index_pix]	= integ[i*3 + j] ;
			}
			free(ipix_final);
			free(localtriad);
			free(integ);
			// this means we will return the TQUmap at the lower resolution nside_variable
			npy_intp npy_shape[2] = {3,npix_variable};
			PyObject *arr 		= PyArray_SimpleNewFromData(2,npy_shape, NPY_DOUBLE, TQUmap);
			PyArray_ENABLEFLAGS((PyArrayObject *)arr, NPY_OWNDATA);
			return(arr);
		}
		else{
			double* TQUmap	= calloc(3*npix_fixed,sizeof(double));
			// else, we are degrading or copying the same
			for (i=0;i<nipix;i++){
				unsigned long index_pix = ipix_final[i];
				// Now I need to assign this pixel to the nside_fixed pixelization
				if (nside_fixed == nside_variable){
					// This means they are at the same pixelization, so it maps directly into TQUmap
					for (j=0;j<3;j++) TQUmap[j*npix_fixed + index_pix]	= integ[i*3 + j] ;
				}
				else if (nside_fixed < nside_variable){
					// This means the filament is smaller than the fixed resolution and we have to degrade the map
					// step will determine how many steps in nside I went, e.g. from 2048 to 8192 I went 2 steps
					unsigned int step = (int)(log(nside_variable)/log(2.0) - log(nside_fixed)/log(2.0)) ;
					unsigned long index_parent_pixel = index_pix >> 2*step ;
					for (j=0;j<3;j++) TQUmap[j*npix_fixed + index_parent_pixel] += integ[i*3 + j]/pow(4,step) ;
					// the division by 4^step is because there are 4^step children pixels inside the parent pixel. We're taking the average. 
				}
			}
			free(ipix_final);
			free(localtriad);
			free(integ);
			// this means we will return the TQUmap at the fixed resolution nside_fixed
			npy_intp npy_shape[2] = {3,npix_fixed};
			PyObject *arr 		= PyArray_SimpleNewFromData(2,npy_shape, NPY_DOUBLE, TQUmap);
			PyArray_ENABLEFLAGS((PyArrayObject *)arr, NPY_OWNDATA);
			return(arr);
		}
	} // this is the end of (isInside==1)
	else{
		printf("Filament %i is outside the box, skipping\n",n_fil) ;
		if (nside_fixed > nside_variable){
			// this means we will return the TQUmap at the lower resolution nside_variable
			double* TQUmap	= calloc(3*npix_variable,sizeof(double));
			npy_intp npy_shape[2] = {3,npix_variable};
			PyObject *arr 		= PyArray_SimpleNewFromData(2,npy_shape, NPY_DOUBLE, TQUmap);
			PyArray_ENABLEFLAGS((PyArrayObject *)arr, NPY_OWNDATA);
			return(arr);
		}
		else{
			// this means we will return the TQUmap at the fixed resolution nside_fixed
			double* TQUmap	= calloc(3*npix_fixed,sizeof(double));
			npy_intp npy_shape[2] = {3,npix_fixed};
			PyObject *arr 		= PyArray_SimpleNewFromData(2,npy_shape, NPY_DOUBLE, TQUmap);
			PyArray_ENABLEFLAGS((PyArrayObject *)arr, NPY_OWNDATA);
			return(arr);
		}
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
