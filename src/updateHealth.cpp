#include "updateHealth.h"

#include <list>

using std::list;

void updateHealth( map<int, vector<double>> &health_parm,
                   const int *ptr_seg, const int &len, const double *ptr_m,
                   const double *ptr_nu2, const double *ptr_intst,
                   const double *ptr_lambda2, const int *ptr_nidx,
                   const double *ptr_nintst, const double *ptr_alpha,
                   const double *ptr_beta, const int &maxit ) {
  // update parm for healthy regions
  for( int i = - 1; i > - 4;  -- i ) {
    int curr_label = i;
    // Rprintf( "curr_label = %d \n", curr_label );
    double mu = -1, sigma2 = 1;
    vector<double> theta;
    for( int  j = 0; j < 6; ++ j ) {
      theta.push_back( 0 );
    }
    list<int> region;
    for( int k = 0; k < len; ++ k ) {
      if( ptr_seg[ 2 * k ] == curr_label ) {
        region.push_back( k + 1 ); // region starts from 1
      }
    }
    int h_idx =  - 1 - curr_label; // == 0, 1, 2
    updateParm( mu, theta, sigma2, region, ptr_m[ h_idx ], 
                ptr_nu2[ h_idx ], ptr_intst, curr_label, 
                ptr_lambda2[ h_idx ], ptr_seg, ptr_nidx,
                ptr_nintst, ptr_alpha[ h_idx ], ptr_beta[ h_idx ],
                                                        maxit );
    
    vector<double> h_parm;
    h_parm.push_back( mu );
    h_parm.push_back( sigma2 );
    
    h_parm.insert( h_parm.end(), theta.begin(), theta.end() );
    health_parm[ curr_label ] = h_parm;
  }
}