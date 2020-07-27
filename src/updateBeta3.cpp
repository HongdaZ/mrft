#include "updateBeta3.h" 
// update beta for t1ce or flair images
void updateBeta3( double *ptr_beta, const double *ptr_alpha, 
                 map<int, vector<double>> &health_parm ) {
  double sigma2;
  for( int i = 0; i < 3; ++ i ){
    sigma2 = health_parm[ - i - 1 ][ 1 ];
    ptr_beta[ i ] = ( ptr_alpha[ i ] + 1 ) * sigma2;
  }
}