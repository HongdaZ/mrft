#include <R.h>
#include <Rinternals.h>

#include "initParm.h"
#include "updateParm.h"
#include "getRegion.h"
#include "clearVector.h"

// Initialize parameters
void initParm( vector<int> &region, vector<double> &theta, 
               const bool first_run, vector<double> &health_parm,
               vector<double> &tumor_parm,
               int *ptr_seg, const double *ptr_m,
               const double *ptr_nu2, const double *ptr_intst, 
               const double *ptr_lambda2, const int *ptr_nidx,
               const double *ptr_nintst, const double *ptr_alpha,
               const double *ptr_beta,
               const vector<int> &n_voxel, const int &n_tumor,
               const double *ptr_a, const double *ptr_b, const int &len, 
               int maxit ) {
  // update parameters for tumor regions
  int count = 0;
  for( int i = 0; i < len; ++ i ) {
    if( n_voxel[ i ] != 0 ) {
      if( first_run || n_voxel[ i ] != 1 ) {
        int curr_label = - i - 4;
        // region starts from 1
        getRegion( region, curr_label, ptr_seg, len );
        double mu = -1, sigma2 = 1; // sigma2 has to be non-zero;
        updateParm( mu, theta, sigma2, region, ptr_m[ 3 ], 
                    ptr_m[ 2 ], ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst,
                    curr_label, ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, 
                    ptr_nintst, ptr_alpha[ 3 ], ptr_beta[ 3 ], maxit );
        int cidx = - curr_label - 4;
        tumor_parm[ 8 * cidx + 0 ] = mu;
        tumor_parm[ 8 * cidx + 1 ] = sigma2;
        for( int j = 0; j < 6; ++ j ) {
          tumor_parm[ 8 * cidx + j + 2 ] = theta[ j ];
        }
      }
      ++ count;
      if( count == n_tumor ) {
        break;
      }
    }
  }
  // Initialize parameters for healthy regions
  for( int i = - 1; i > - 4;  -- i ) {
    int curr_label = i;
    // Rprintf( "curr_label = %d \n", curr_label );
    double mu = -1, sigma2 = 1;
    getRegion( region, curr_label, ptr_seg, len );
    int h_idx =  - 1 - curr_label; // == 0, 1, 2
    updateParm( mu, theta, sigma2, region, ptr_m[ h_idx ], 
                ptr_nu2[ h_idx ], ptr_intst, curr_label, 
                ptr_lambda2[ h_idx ], ptr_seg, ptr_nidx, ptr_nintst,
                ptr_alpha[ h_idx ], ptr_beta[ h_idx ], maxit );
    
    int cidx = - curr_label - 1;
    health_parm[ 8 * cidx + 0 ] = mu;
    health_parm[ 8 * cidx + 1 ] = sigma2;
    for( int j = 0; j < 6; ++ j ) {
      health_parm[ 8 * cidx + j + 2 ] = theta[ j ];
    }
  }
  return;
}
// Initialize parameters for t1ce and flair images
void initParmHealth3( vector<int> &region, vector<double> &theta,
                      vector<double> &health_parm,
                      int *ptr_seg, const double *ptr_m,
                      const double *ptr_nu2, const double *ptr_intst, 
                      const double *ptr_lambda2, const int *ptr_nidx,
                      const double *ptr_nintst, const double *ptr_alpha,
                      const double *ptr_beta, const int &len, 
                      int maxit ) {
  // Initialize parameters for healthy regions
  for( int i = - 1; i > - 4;  -- i ) {
    int curr_label = i;
    double mu = -1, sigma2 = 1;
    getRegion( region, curr_label, ptr_seg, len );
    int h_idx =  - 1 - curr_label; // == 0, 1, 2
    updateParm( mu, theta, sigma2, region, ptr_m[ h_idx ], 
                ptr_nu2[ h_idx ], ptr_intst, curr_label,
                ptr_lambda2[ h_idx ], ptr_seg, ptr_nidx, ptr_nintst, 
                ptr_alpha[ h_idx ], ptr_beta[ h_idx ], maxit );
    
    int cidx = - curr_label - 1;
    health_parm[ 8 * cidx + 0 ] = mu;
    health_parm[ 8 * cidx + 1 ] = sigma2;
    for( int j = 0; j < 6; ++ j ) {
      health_parm[ 8 * cidx + j + 2 ] = theta[ j ];
    }
  }
  return;
}
