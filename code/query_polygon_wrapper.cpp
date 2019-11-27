#include <healpix_base.h>
#include <rangeset.h>
#include <pointing.h>
#include <vec3.h>
#include <vector>
#include <cmath>

extern "C" {
	void query_polygon_wrapper(double theta[4], double phi[4], long nside, long* ipix_arr, long n_ipix){
		// theta and phi are array with size 4
		pointing a,b,c,d;
		a.theta=theta[0]; b.theta=theta[1];c.theta=theta[2];d.theta=theta[3];
		a.phi=phi[0];b.phi=phi[1];c.phi=phi[2];d.phi=phi[3];
		std::vector<pointing> vertices;

		vertices.push_back(a);
		vertices.push_back(b);
		vertices.push_back(c);
		vertices.push_back(d);
		
		T_Healpix_Base<long> hp_base(nside,RING,SET_NSIDE); 
		rangeset<long> ipix;	
		hp_base.query_polygon(vertices,ipix);
		std::vector<long> v = ipix.toVector();
		
		n_ipix	= v.size();

		for(std::size_t i = 0; i < v.size(); ++i) {
			//std::cout << v[i] << "\n";
			ipix_arr[i]	= v[i];
		}
		return;
	}
}
