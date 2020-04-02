double** FilamentPaint_RotationMatrix(double angles_arr[2]);
double** FilamentPaint_InvertRotMat(double rot_matrix[3][3]);
double** FilamentPaint_xyzVertices(double rot_matrix[3][3], double sizes_arr[3], double centers_arr[3], double Size, int* isInside);
double** FilamentPaint_xyzNormalToFaces(double rot_matrix[3][3]);
double*** FilamentPaint_xyzFaces(double xyz_vertices[8][3]);
double*** FilamentPaint_xyzEdgeVectors(double xyz_faces[6][4][3]);
double*** FilamentPaint_xyzEdgeVectorsUnit(double xyz_edges[6][2][3]);
void FilamentPaint_DoQueryPolygon(long nside_long,double xyz_faces[6][4][3], long* ipix, long* nipix);
double* FilamentPaint_CalculateDistances(double xyz_normal_to_faces[6][3], double xyz_faces[6][4][3], double xyz_edges[6][2][3], double xyz_edges_unit[6][2][3], double rUnitVector_ipix[3]);
double FilamentPaint_Density(double r, double rot_matrix[3][3], double rUnitVector_ipix[3], double centers_arr[3], double sizes_arr[3]);
double* FilamentPaint_TrilinearInterpolation(PyObject* Bcube_obj, double size_box, int nbox, double vector[3]);
double* FilamentPaint_Bxyz(double r, PyObject* Bcube_obj,double size_box, int nbox,double rUnitVector_ipix[3],double local_triad[3][3]);
double* FilamentPaint_Integrator(double r1, double r2, double rot_matrix[3][3], double rUnitVector_ipix[3],double local_triad[3][3], double centers_arr[3], double sizes_arr[3], PyObject* Bcube_obj, double size_box, int nbox);
double* FilamentPaint_Blocal(PyObject* Bcube_obj,double size_box, int nbox,double center_arr[3]);