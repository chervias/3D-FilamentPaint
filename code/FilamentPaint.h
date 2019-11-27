double** FilamentPaint_RotationMatrix(double angles_arr[2]);
double** FilamentPaint_xyzVertices(double rot_matrix[3][3], double sizes_arr[3], double centers_arr[3]);
double** FilamentPaint_xyzNormalToFaces(double rot_matrix[3][3]);
double*** FilamentPaint_xyzFaces(double xyz_vertices[8][3]);
double*** FilamentPaint_xyzEdgeVectors(double xyz_faces[6][4][3]);
