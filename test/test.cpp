// Header file for input output functions 
#include <healpix_base.h>
#include <rangeset.h>
#include <pointing.h>
#include <vec3.h>
#include <vector>
  
using namespace std; 
  
int main(){ 
	// prints hello world 
	//cout<<"Hello World\n";

	//vec3_t<float64> b 	= vec3_t<float64>(-0.1,0.1,1.0) ;
	//vec3_t<float64> c 	= vec3_t<float64>(-0.1,-0.1,1.0) ;
	//vec3_t<float64> d 	= vec3_t<float64>(0.1,-0.1,1.0) ;
	
	pointing aa;
	aa.theta 	=0.0;
	aa.phi		=0.0;
	
	std::vector<pointing> vertices;
	vertices.push_back(aa);
	vertices.push_back(aa);
	vertices.push_back(aa);
	vertices.push_back(aa);
	
	std::cout << "Size: " << vertices.size() << '\n';
	
	

	return 0; 
} 
