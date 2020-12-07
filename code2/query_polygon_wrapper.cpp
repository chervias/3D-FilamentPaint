#include <healpix_base.h>
#include <rangeset.h>
#include <pointing.h>
#include <vec3.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

extern "C" {
	void query_polygon_wrapper(double* theta, double* phi, unsigned int nside, unsigned long* ipix_arr, unsigned long* nipix){
		// theta and phi are array with size 4
		pointing a,b,c,d;
		a.theta=theta[0]; b.theta=theta[1];c.theta=theta[2];d.theta=theta[3];
		a.phi=phi[0];b.phi=phi[1];c.phi=phi[2];d.phi=phi[3];
		std::vector<pointing> vertices;
		
		vertices.push_back(a);
		vertices.push_back(b);
		vertices.push_back(c);
		vertices.push_back(d);
		
		//printf("nside = %i \n",nside);
		//printf("first = (%.20E , %.20E) \n",a.theta,a.phi);
		//printf("second = (%.20E ,%.20E) \n",b.theta,b.phi);
		//printf("third = (%.20E ,%.20E) \n",c.theta,c.phi);
		//printf("fourth = (%.20E , %.20E) \n",d.theta,d.phi);
		
		T_Healpix_Base<long> hp_base(nside,NEST,SET_NSIDE); 
		rangeset<long> ipix;	
		hp_base.query_polygon(vertices,ipix);
		std::vector<long> v = ipix.toVector();
		*nipix	= v.size();
		for(std::size_t i = 0; i < v.size(); i++) {
			ipix_arr[i]	= v[i];
		}
	}
}

extern "C" {
	void localtriad_wrapper(unsigned long* ipix_final, unsigned long nipix, unsigned int nside, double* localtriad){
		// localtriad will be nipix by 3 (hat(x),hat(y),hat(z)) by 3 (x,y,z components) matrix
		// localtriad[nipix * 3 * 3] = local_triad[i][j][k] = localtriad[i*3*3 + j*3 + k ]
		// set the healpix base object
		T_Healpix_Base<long> hp_base(nside,NEST,SET_NSIDE);
		for(std::size_t i = 0; i < nipix; i++){
			vec3 vec = hp_base.pix2vec(ipix_final[i]);
			double theta = atan2(sqrt(pow(vec.x,2)+pow(vec.y,2)),vec.z);
			double phi = atan2(vec.y,vec.x);
			double norm = sqrt(pow(vec.x,2)+pow(vec.y,2)+pow(vec.z,2)) ;
			//std::cout << "Pixel id=" << ipix[i] << " Vector |" << vec.x << "|" << vec.y << "|" << vec.z << " Angles " << theta << phi << "\n";
			// fill the hat(x) vector from local triad (equivalent to hat(theta) from spherical coordinates)
			localtriad[i*3*3 + 0*3 + 0] = cos(theta)*cos(phi);
			localtriad[i*3*3 + 0*3 + 1] = cos(theta)*sin(phi);
			localtriad[i*3*3 + 0*3 + 2] = -sin(theta);
			// fill the hat(y) vector from local triad (equivalent to hat(phi) from spherical coordinates)
			localtriad[i*3*3 + 1*3 + 0] = -sin(phi);
			localtriad[i*3*3 + 1*3 + 1] = cos(phi);
			localtriad[i*3*3 + 1*3 + 2] = 0.0;
			// fill the hat(z) vector from local triad (equivalent to hat(r) from spherical coordinates or the output from pix2vec)
			localtriad[i*3*3 + 2*3 + 0] = vec.x/norm;
			localtriad[i*3*3 + 2*3 + 1] = vec.y/norm;
			localtriad[i*3*3 + 2*3 + 2] = vec.z/norm;
		}
	}
}

/*
extern "C" {
	void upgrade_with_kernel(double fwhm, unsigned int nside_low, unsigned int nside_high, unsigned int npix_high, unsigned int nipix, unsigned int* ipix_final,double* integ, double* TQUmap){
		// this function will upgrade a filament sampled at low resolution to a higher resolution using a 2d gaussian kernel
		// inputs: 
		//fwhm of the kernel
		//nside low reso
		//nside high reso
		//nipix is the number of pixels in the low resolution map
		// ipix_final --> the indices of the pixels that are within the filament.
		//integ --> the values of tqu for each of the pixel in ipix_final
		// outputs that will be modified: 
		// TQUmap 
		// Code
		double sigma = fwhm/2.355 ;
		unsigned int step = (int)(log(nside_high)/log(2.0) - log(nside_low)/log(2.0)) ;
		T_Healpix_Base<long> hp_base_low(nside_low,NEST,SET_NSIDE);
		T_Healpix_Base<long> hp_base_high(nside_high,NEST,SET_NSIDE);
		// We iterate over the low reso pixels and store the angles of the centers so we do it only once
		//double* theta_low = static_cast<double*>(calloc(nipix,sizeof(double)));
		//double* phi_low = static_cast<double*>(calloc(nipix,sizeof(double)));
		//for(std::size_t i = 0; i < nipix; i++){
		//	std::cout << "Low reso " << i << "/" << ipix_final[i] << "\n";
		//	vec3 vec     = hp_base_low.pix2vec(ipix_final[i]);
		//	theta_low[i] = atan2(sqrt(pow(vec.x,2)+pow(vec.y,2)),vec.z);
		//	phi_low[i]   = atan2(vec.y,vec.x);		
		//}
		// We have to iterate over the children pixels in the high resolution map
		// for that we iterate over the low resolution pixels and then over the children pixels
		for(std::size_t i = 0; i < nipix; i++){
			unsigned int index_pix = ipix_final[i];
			unsigned int index_children_pixel = index_pix << 2*step ;
			for(std::size_t j=index_children_pixel; j<(index_children_pixel+pow(4,step)); j++){
				// here we are iterating over the children pixels
				// we need the theta,phi of the children pixel
				vec3 vec = hp_base_high.pix2vec(j);
				double theta_high = atan2(sqrt(pow(vec.x,2)+pow(vec.y,2)),vec.z);
				double phi_high   = atan2(vec.y,vec.x);
				// Now we need to know the low reso pixels that belong within 5 sigma
				pointing ctr;
				ctr.theta = theta_high;
				ctr.phi = phi_high;
				rangeset<long> ipix; 
				hp_base_low.query_disc(ctr,10*sigma,ipix);
				// now ipix contains the low reso pixels over which I will iterate to sum
				double norm = 0.0 ;
				std::vector<long> v = ipix.toVector();
				//std::cout << "Pixel " << j << "has " << v.size() << "big pixels around it \n" ;
				// Now we iterate over the low reso pixels, where we will sum with a weight given by the gaussian kernel
				std::vector<int> myvector (ipix_final,ipix_final+nipix);
				for(std::size_t k = 0; k < v.size(); k++){
					// I need to find the index of v[k] in the ipix_final array (and therefore on the integ array)
					unsigned int elem = v[k];
					//std::cout << elem << "\n" ;
					std::vector<int>::iterator pp;
					pp = find (myvector.begin(), myvector.end(), elem);
					vec3 vec     = hp_base_low.pix2vec(elem);
					double theta_low = atan2(sqrt(pow(vec.x,2)+pow(vec.y,2)),vec.z);
					double phi_low   = atan2(vec.y,vec.x);
					// distance between the high reso pixel and the low reso pixel, with haversines formula
					// theta = 90 deg - lat
					// ang 1 will be low reso , ang 2 will be high
					double distance = 2*asin(sqrt( pow(sin(0.5*(-theta_high+theta_low)),2) + cos(M_PI/2 - theta_high)*cos(M_PI/2 - theta_low)*pow(sin(0.5*(phi_high-phi_low)),2) )) ;
					double weight = exp(-0.5*pow(distance,2)/pow(sigma,2));
					norm = norm + weight ;
					// now if the pixel elem was not found on ipix_final, it means it was in the edges of the filament, so its value on integ is 0
					if (pp != myvector.cend()){
						// this means elem was found in ipix_final
						//std::cout << "Pixel " << elem << " " << ipix_final[std::distance(myvector.begin(), pp)] << " in position " << std::distance(myvector.begin(), pp) << "\n";
						// now we sum the corresponding low reso pixel to the high reso pixel 
						// index l loops over t,q,u
						unsigned int index_ = std::distance(myvector.begin(), pp);
						for (std::size_t l=0;l<3;l++) TQUmap[l*npix_high + j] += weight*integ[index_*3+l] ;
					}
					//else
						// this is the pixel elem  is not present in ipix_final, therefore the value in integ is 0
						// we sum 0 to the corresponding pixel in TQUmap, i.e. we do nothing
				}
				//std::cout << "Pixel " << j << "has normalization " << norm << "\n" ;
				for (std::size_t l=0;l<3;l++) TQUmap[l*npix_high + j] =  TQUmap[l*npix_high + j] / norm ;
			}
		}
	}
}
*/