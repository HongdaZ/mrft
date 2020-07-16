#include "updateBeta.h"
// update beta
void updateBeta( double *ptr_beta, const double *ptr_alpha, 
                 map<int, vector<double>> &health_parm,
                 map<int, list<int>> &tumor_regions,
                 map<int, vector<double>> &tumor_parm ) {
  
  for( int i = 0; i < 3; ++ i ){
    double sigma2 = health_parm[ - i - 1 ][ 1 ];
    ptr_beta[ i ] = ( ptr_alpha[ i ] + 1 ) * sigma2;
  }
  // update beta for tumor regions based on the date of
  // the biggest tumor region
  int max_size = 0;
  int max_label = 0;
  for( map<int, list<int>>::iterator it = tumor_regions.begin();
       it != tumor_regions.end(); ++ it ) {
    if( it->second.size() > max_size ) {
      max_size = it->second.size();
      max_label = it->first;
    }
  }
  double tumor_sigma2 = tumor_parm[ max_label ][ 1 ];
  ptr_beta[ 3 ] = ( ptr_alpha[ 3 ] + 1 ) * tumor_sigma2;
}