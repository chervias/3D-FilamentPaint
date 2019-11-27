// Header file for input output functions 
#include <healpix_base.h>
#include <rangeset.h>
#include <pointing.h>
#include <vec3.h>
#include <vector>
#include <cmath>

using namespace std; 
  
int main(){ 	
	pointing a,b,c,d;
	a.theta 	= 0.5*M_PI - 0.01;
	a.phi		= 0.01;
	
	b.theta 	= 0.5*M_PI - 0.01;
	b.phi		= -0.01;
	
	c.theta 	= 0.5*M_PI + 0.01;
	c.phi		=-0.1;
	
	d.theta 	= 0.5*M_PI + 0.01;
	d.phi		=+0.1;
	
	vector<pointing> vertices;
	
	vertices.push_back(a);
	vertices.push_back(b);
	vertices.push_back(c);
	vertices.push_back(d);
	
	int nside=512;
	
	T_Healpix_Base<int> hp_base(nside,RING,SET_NSIDE);
	rangeset<int> ipix;		
	hp_base.query_polygon(vertices,ipix);
	
	vector<int> v = ipix.toVector();
	
	// Using a for loop with index
	for(std::size_t i = 0; i < v.size(); ++i) {
		std::cout << v[i] << "\n";
	}
	
	//std::cout << "pixels: " << pixels << '\n';
	
	

	return 0; 
} 
