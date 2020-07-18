#include <R.h>
#include <Rinternals.h>

#include <cmath>   

#include "updateParm.h"
#include "labelRegion.h" 

using std::abs;

// previous updateParm is slow because of using std::find with std::list
// and extra times of copying matrices

// update parameters for healthy cells
void updateParm( double &mu, vector<double> &theta, double &sigma2, 
                 const list<int> &region,
                 const double m,
                 const double nu2,
                 const double *ptr_intst,
                 int curr_label,
                 const double lambda2,
                 int *ptr_seg,
                 const int *ptr_nidx,
                 const double *ptr_nintst,
                 const double alphal,
                 const double betal,
                 int maxit ) {

  for( int i = 0; i < 6; ++ i ) {
    theta[ i ] = 0;
  }
  int i = 0;
  double tol = 1;
  double tmp = 0;
  mu = - 1;
  // Initialize matrices and vectors
  int nrow = region.size();
  int ncol = 6;
  double *yln_ = new double[ nrow * ncol ];
  double *yln_i = new double[ nrow * ncol ];
  double *yl_ = new double[ nrow ];
  
  list<int>::const_iterator it = region.begin();
  int idx;
  double sum_y = 0;
  for( int j = 0; j < nrow; ++ it, ++ j ) {
    idx = *it;
    yl_[ j ] = ptr_intst[ idx - 1 ];
    sum_y += yl_[ j ];
  }
  for( int i = 0; i < 6; ++ i ) {
    it = region.begin();
    for( int j = 0; j < nrow; ++ it, ++ j ) {
      idx = *it;
      int nidx =  ptr_nidx[ 6 * ( idx - 1 ) + i ];
      if( nidx != NA_INTEGER ) {
        // healthy
        if( curr_label < 0 && curr_label > - 4 ) {
          int nlabel = ptr_seg[ 2 * ( nidx - 1 ) ];
          if( nlabel == curr_label ) {
            yln_[ i * nrow + j ] = ptr_nintst[ 6 * ( idx - 1 ) + i ];
            yln_i[ i * nrow + j ] = 1;
          } else {
            yln_[ i * nrow + j ] = 0;
            yln_i[ i * nrow + j ] = 0;
          }
          // tumor
        } else {
          int nlabel = ptr_seg[ 2 * nidx - 1 ];
          if( nlabel == 3 ) {
            yln_[ i * nrow + j ] = ptr_nintst[ 6 * ( idx - 1 ) + i ];
            yln_i[ i * nrow + j ] = 1;
          } else {
            yln_[ i * nrow + j ] = 0;
            yln_i[ i * nrow + j ] = 0;
          }
        } 
      } else {
        yln_[ i * nrow + j ] = 0;
        yln_i[ i * nrow + j ] = 0;
      }
    }
  }
  double *yln = new double[ nrow * ncol ];
  double *yl = new double[ nrow ];
  
  double sum_theta;
  while( i < 2 || ( i < maxit && tol > .0001 )  ) {
    sum_theta = 0;
    for( int j = 0; j < 6; ++ j ) {
      sum_theta += theta[ j ];
    }
    tmp = updateMu( nrow, sigma2, m, nu2, sum_theta, sum_y );
    // Rprintf( "mu = %f \n", tmp );
    tol = abs( mu - tmp );
    mu = tmp;
    updateTS( nrow, curr_label, mu, sigma2, lambda2, ptr_seg, ptr_nidx,
              ptr_intst, ptr_nintst, theta, alphal, betal, yl, yln,
              yln_, yln_i, yl_ );
    // for( int j = 0; j < 6; ++ j ) {
    //   Rprintf( "%f\t", theta[ j ] );
    // }
    // Rprintf( "\n" );
    // Rprintf( "sigma2 = %f\n", sigma2 );
    ++ i;
  }
  delete [] yln_;
  delete [] yln_i;
  delete [] yl_;
  delete [] yln;
  delete [] yl;
  
  return;
}

// update parameters for tumor cells
void updateParm( double &mu, vector<double> &theta, double &sigma2, 
                 const list<int> &region,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double *ptr_intst,
                 const int curr_label,
                 const double lambda2,
                 int *ptr_seg,
                 const int *ptr_nidx,
                 const double *ptr_nintst,
                 const double alphal,
                 const double betal,
                 int maxit ) {
  for( int i = 0; i < 6; ++ i ) {
    theta[ i ] = 0;
  }
  int i = 0;
  double tol = 1;
  double tmp = 0;
  mu = - 1;
  // the voxel in region is labelled as 3 in ptr_seg[ 2, ]
  labelRegion( region, ptr_seg );
  
  // Initialize matrices and vectors
  int nrow = region.size();
  int ncol = 6;
  double *yln_ = new double[ nrow * ncol ];
  double *yln_i = new double[ nrow * ncol ];
  double *yl_ = new double[ nrow ];
  
  list<int>::const_iterator it = region.begin();
  int idx;
  double sum_y = 0;
  for( int j = 0; j < nrow; ++ it, ++ j ) {
    idx = *it;
    yl_[ j ] = ptr_intst[ idx - 1 ];
    sum_y += yl_[ j ];
  }
  for( int i = 0; i < 6; ++ i ) {
    it = region.begin();
    for( int j = 0; j < nrow; ++ it, ++ j ) {
      idx = *it;
      int nidx =  ptr_nidx[ 6 * ( idx - 1 ) + i ];
      if( nidx != NA_INTEGER ) {
        // healthy
        if( curr_label < 0 && curr_label > - 4 ) {
          int nlabel = ptr_seg[ 2 * ( nidx - 1 ) ];
          if( nlabel == curr_label ) {
            yln_[ i * nrow + j ] = ptr_nintst[ 6 * ( idx - 1 ) + i ];
            yln_i[ i * nrow + j ] = 1;
          } else {
            yln_[ i * nrow + j ] = 0;
            yln_i[ i * nrow + j ] = 0;
          }
          // tumor
        } else {
          int nlabel = ptr_seg[ 2 * nidx - 1 ];
          if( nlabel == 3 ) {
            yln_[ i * nrow + j ] = ptr_nintst[ 6 * ( idx - 1 ) + i ];
            yln_i[ i * nrow + j ] = 1;
          } else {
            yln_[ i * nrow + j ] = 0;
            yln_i[ i * nrow + j ] = 0;
          }
        } 
      } else {
        yln_[ i * nrow + j ] = 0;
        yln_i[ i * nrow + j ] = 0;
      }
    }
  }
  double *yln = new double[ nrow * ncol ];
  double *yl = new double[ nrow ];

  double sum_theta;
  while(  i < 2 || ( i < maxit && tol > .0001 ) ) {
    sum_theta = 0;
    for( int j = 0; j < 6; ++ j ) {
      sum_theta += theta[ j ];
    }
    tmp = updateMu( nrow, sigma2, m, mk_1, a, b, sum_theta, sum_y );
    // Rprintf( "mu = %f \n", tmp );
    tol = abs( mu - tmp );
    mu = tmp;
    updateTS( nrow, curr_label, mu, sigma2, lambda2, ptr_seg, ptr_nidx,
              ptr_intst, ptr_nintst, theta, alphal, betal, yl, yln,
              yln_, yln_i, yl_ );
    // for( int j = 0; j < 6; ++ j ) {
    //   Rprintf( "%f\t", theta[ j ] );
    // }
    // Rprintf( "\n" );
    // Rprintf( "sigma2 = %f\n", sigma2 );
    ++ i;
  }
  // change ptr_seg[ 2, region ] back
  recoverLabel( region, ptr_seg );
  delete [] yln_;
  delete [] yln_i;
  delete [] yl_;
  delete [] yln;
  delete [] yl;
  return;
}

// update parameters for outliers
void updateParm( double &mu, double &sigma2, 
                 const list<int> &region,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double *ptr_intst,
                 int curr_label,
                 const double lambda2,
                 int *ptr_seg,
                 const int *ptr_nidx,
                 const double alphal,
                 const double betal,
                 int maxit ) {
  int i = 0;
  double tol = 1;
  double tmp = 0;
  mu = - 1;
  list<int>::const_iterator it = region.begin();
  int idx = *it;
  double sum_y = ptr_intst[ idx - 1 ];
  int n = 1;
  while(  i < 2 || ( i < maxit && tol > .0001 ) ) {
    tmp = updateMu( n, sigma2, m, mk_1, a, b, sum_y );
    // Rprintf( "mu = %f \n", tmp );
    tol = abs( mu - tmp );
    mu = tmp;
    updateSigma( idx, mu, ptr_intst, alphal, betal, sigma2 );
    // Rprintf( "\n" );
    // Rprintf( "sigma2 = %f\n", sigma2 );
    ++ i;
  }
  return;
}