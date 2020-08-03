#include "cmpE.h"
#include "label2col.h"
#include "getParm.h"

// compare energy for t1ce or flair images in estimation
void cmpE3( int idx, vector<double> &health_parm,
            int *ptr_seg, const int *ptr_nidx,
            const double *ptr_intst, const double *ptr_nintst,
            const double *ptr_delta, const double *ptr_gamma,
            vector<double> &theta ) {
  
  double energy;
  double min_energy;
  int min_label;
  double mu;
  double sigma2;
  int idxcol;
  for( int i = - 1; i > - 4; -- i ) {
    idxcol = label2col( i );
    getParm( mu, sigma2, theta, health_parm, idxcol );
    energy = energyY( i, idx, mu, sigma2, ptr_seg, ptr_nidx,
                      ptr_intst, ptr_nintst, theta );
    energy += energyX( i, idx, false, ptr_seg, ptr_nidx,
                       ptr_delta[ 0 ], ptr_gamma[ 0 ] );
    if( i == - 1 ) {
      min_energy = energy;
      min_label = i;
    } else {
      if( energy < min_energy ) {
        min_energy = energy;
        min_label = i;
      }
    }
  }
  ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
}