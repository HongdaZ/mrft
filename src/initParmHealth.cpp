#include <R.h>
#include <Rinternals.h>

#include "initParmHealth.h"
#include "updateParm.h"

// Initialize parameters for t1ce and flair images
void initParmHealth3( map<int, vector<double>> &health_parm,
                     int *ptr_seg, const double *ptr_m,
                     const double *ptr_nu2, const double *ptr_intst, 
                     const double *ptr_lambda2, const int *ptr_nidx,
                     const double *ptr_nintst, const double *ptr_alpha,
                     const double *ptr_beta, int len, 
                     int maxit ) {
  // Initialize parameters for healthy regions
  for( int i = - 1; i > - 4;  -- i ) {
    int curr_label = i;
    double mu = -1, sigma2 = 1;
    vector<double> theta( 6, 0 );
    list<int> region;
    for( int k = 0; k < len; ++ k ) {
      if( ptr_seg[ 2 * k ] == curr_label ) {
        region.push_back( k + 1 ); // region starts from 1
      }
    }
    int h_idx =  - 1 - curr_label; // == 0, 1, 2
    updateParm( mu, theta, sigma2, region, ptr_m[ h_idx ], ptr_nu2[ h_idx ],
                ptr_intst, curr_label, ptr_lambda2[ h_idx ], ptr_seg, ptr_nidx,
                ptr_nintst, ptr_alpha[ h_idx ], ptr_beta[ h_idx ], maxit );
    
    vector<double> h_parm( 2, 0 );
    h_parm[ 0 ] = mu;
    h_parm[ 1 ] = sigma2;
    
    h_parm.insert( h_parm.end(), theta.begin(), theta.end() );
    health_parm[ curr_label ] = h_parm;
  }
  return;
}
