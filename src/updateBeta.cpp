#include "updateBeta.h"
// update beta for healthy and tumor regions
void updateBeta( double *ptr_beta, const double *ptr_alpha, 
                 const vector<double> &health_parm,  
                 const vector<int> &n_voxel, const int &n_tumor,
                 const vector<double> &tumor_parm ) {
  double sigma2;
  for( int i = 0; i < 3; ++ i ) {
    sigma2 = health_parm[ i * 8 + 1 ];
    ptr_beta[ i ] = ( ptr_alpha[ i ] + 1 ) * sigma2;
  }
  // update beta for tumor regions based on the date of
  // the biggest tumor region
  int max_size = 0;
  int max_label = 0;
  int count = 0;
  int len = n_voxel.size();
  for( int i = 0; i < len; ++ i ) {
    if( n_voxel[ i ] > 0 ) {
      if( n_voxel[ i ] > max_size ) {
        max_size = n_voxel[ i ];
        max_label = - i - 4;
      }
      ++ count;
      if( count == n_tumor ) {
        break;
      }
    }
  }
  int cidx = - max_label - 4;
  double tumor_sigma2 = tumor_parm[ 8 * cidx + 1 ];
  ptr_beta[ 3 ] = ( ptr_alpha[ 3 ] + 1 ) * tumor_sigma2;
}
// update beta for t1ce or flair images
void updateBeta3( double *ptr_beta, const double *ptr_alpha, 
                  const vector<double> &health_parm ) {
  double sigma2;
  for( int i = 0; i < 3; ++ i ) {
    sigma2 = health_parm[ i * 8 + 1 ];
    ptr_beta[ i ] = ( ptr_alpha[ i ] + 1 ) * sigma2;
  }
}