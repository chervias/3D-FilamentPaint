#include "healpix_map.h"

extern "C" {

  void hp_interp_ring_cxx_(double *theta, double *phi, int *npt, double *themap, int *nside, double *val) {

  int npix = 12* (*nside) * (*nside);
  
 
  Healpix_Map<double> hp_map((*nside),RING,SET_NSIDE);
  int i;

  for (i=0;i<npix;i++) {
    hp_map[i] = themap[i];
  }
  
  pointing pt;
  for (i=0;i< *npt;i++) {
    pt.theta = theta[i];
    pt.phi = phi[i];
    val[i] =  hp_map.interpolated_value2(pt);
  }
 
  return;
}


}
