#include <R.h>
#include <Rinternals.h>

#include "initParm.h"
#include "updateParm.h"
#include "energyY.h"
#include "energyX.h"

// Initialize parameters
void initParm( map<int, vector<double>> &health_parm,
               map<int, vector<double>> &tumor_parm,
               const int *ptr_seg, const double *ptr_m,
               const double *ptr_nu2, const double *ptr_intst, 
               const double *ptr_lambda2, const int *ptr_nidx,
               const double *ptr_nintst, const double *ptr_alpha,
               const double *ptr_beta,
               map<int, list<int>> &tumor_regions, 
               const double *ptr_a, const double *ptr_b, int len, 
               int maxit ) {
  // update parameters for tumor regions
  for( map<int, list<int>>::iterator it = tumor_regions.begin();
       it != tumor_regions.end(); ++ it ) {
    int curr_label = it->first;
    list<int> &region = it->second;
    double mu = -1, sigma2 = 1; // sigma2 has to be non-zero;
    vector<double> theta( 6, 0 );
    updateParm( mu, theta, sigma2, region, ptr_m[ 3 ], ptr_m[ 2 ],
                ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst, curr_label, 
                ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, ptr_nintst,
                ptr_alpha[ 3 ], ptr_beta[ 3 ], maxit );
    vector<double> t_parm( 8, 0 );
    t_parm[ 0 ] = mu;
    t_parm[ 1 ] = sigma2;
    for( int i = 0; i < 6; ++ i ) {
      t_parm[ i + 2 ] = theta[ i ];
    }
    tumor_parm[ curr_label ] = t_parm;
    Rprintf( "label = %d; mu = %f; sigma2 = %f; theta = ", curr_label, mu,
             sigma2 );
    for( int i = 0; i < 6; ++ i ) {
      Rprintf( "%f, ", theta[ i ] );
    }
    Rprintf( "\n" );
  }
  // Initialize parameters for healthy regions
  for( int i = - 1; i > - 4;  -- i ) {
    int curr_label = i;
    // Rprintf( "curr_label = %d \n", curr_label );
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
    Rprintf( "label = %d; mu = %f; sigma2 = %f; theta = ", curr_label, mu,
             sigma2 );
    for( int i = 0; i < 6; ++ i ) {
      Rprintf( "%f, ", theta[ i ] );
    }
    Rprintf( "\n" );
  }
  return;
}
