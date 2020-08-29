#include "updateBeta.h"
#include "label2col.h"

// update beta for healthy and tumor regions
void updateBeta( double *ptr_beta, const double *ptr_alpha, 
                 const vector<double> &health_parm,  
                 const list<list<int>> &tumor_regions,
                 const vector<double> &tumor_parm ) {
  double sigma2;
  for( int i = 0; i < 3; ++ i ) {
    sigma2 = health_parm[ i * 8 + 1 ];
    ptr_beta[ i ] = ( ptr_alpha[ i ] + 1 ) * sigma2;
  }
  // update beta for tumor regions based on the parameter of
  // the biggest tumor region
  int max_size = 0;
  int max_label = 0;
  int curr_size = 0;
  for( list<list<int>>::const_iterator it = tumor_regions.begin();
       it != tumor_regions.end(); ++ it ) {
    curr_size = it->size() - 1;
    if( curr_size > max_size ) {
      max_size = curr_size;
      max_label = it->front();
    }
  }
  int cidx = label2col( max_label );
  sigma2 = tumor_parm[ 8 * cidx + 1 ];
  ptr_beta[ 3 ] = ( ptr_alpha[ 3 ] + 1 ) * sigma2;
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
// update beta for t2 images
void updateBeta2( double *ptr_beta, const double *ptr_alpha, 
                  const vector<double> &health_parm ) {
  double sigma2;
  for( int i = 0; i < 2; ++ i ) {
    sigma2 = health_parm[ i * 8 + 1 ];
    ptr_beta[ i ] = ( ptr_alpha[ i ] + 1 ) * sigma2;
  }
}