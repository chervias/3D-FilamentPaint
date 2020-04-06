#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
 
int main(){
	int i,j;
	double determinant=0;
	double angles[2] = {0.6,0.5};
	double** a = FilamentPaint_RotationMatrix(angles) ;
 
  for(i=0;i<3;i++)
      determinant = determinant + (a[0][i]*(a[1][(i+1)%3]*a[2][(i+2)%3] - a[1][(i+2)%3]*a[2][(i+1)%3]));
 
   printf("\nInverse of matrix is: \n\n");
   for(i=0;i<3;i++){
      for(j=0;j<3;j++)
           printf("%.2f\t",((a[(i+1)%3][(j+1)%3] * a[(i+2)%3][(j+2)%3]) - (a[(i+1)%3][(j+2)%3]*a[(i+2)%3][(j+1)%3]))/ determinant);
       printf("\n");
   }
 
   return 0;
}