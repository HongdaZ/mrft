#include <R.h>
#include <Rinternals.h>

#include <vector>

#include "energyX.h"

using std::vector;

// pairwise potential function
double pairwise( const int &left, const int &right ) {
  double ans = 1;
  if( left == right || left >= 0 || right >= 0 ) {
    ans = 0;
  }
  return ans;
}

double unary( const int &label ) {
  double ans = 0;
  if( label > 0 ) {
    ans = 1;
  }
  return ans;
}

// region starts from 1
// calculate energyX for single voxel
double energyX( const int &curr_label, const int &curr_idx, 
                const bool &region,
                const int *ptr_seg, const int *ptr_nidx,
                const double &delta, const double &gamma ) {
  double energy = delta * unary( curr_label );
  int nbr_idx = 0;
  int label;
  double scale;
  for( int i = 0; i < 6; ++ i ) {
    scale = 1;
    nbr_idx = ptr_nidx[ 6 * ( curr_idx - 1 ) + i ];
    if( nbr_idx != NA_INTEGER ) {
      label = ptr_seg[ 2 * ( nbr_idx - 1 ) ];
      if( region && label < - 3 ) {
        scale = 2;
      }
      energy += scale * gamma * pairwise( curr_label, label );
    }
  }
  return energy;
}
