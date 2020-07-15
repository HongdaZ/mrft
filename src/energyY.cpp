#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include <set>
#include <vector>

#include "energyY.h"

using std::set;
using std::vector;



// region starts from 1
// calculate energy for tumor regions
double energyY( const set<int> &region,
               double mu,
               double mk1,
               double sigma2,
               double lambda2,
               const int *ptr_seg,
               const int *ptr_nidx,
               const double *ptr_intst,
               const double *ptr_nintst,
               const vector<double> &theta,
               double alphak,
               double betak,
               double a, double b ){
  int nrow = region.size();
  int ncol = 6;
  double *yln = new double[ nrow * ncol ];
  double *yl = new double[ nrow ];
  
  set<int>::iterator it = region.begin();
  int idx;
  // Initialize yln and yl;
  for( int j = 0; it!= region.end(); ++ it, ++ j ) {
    idx = *it;
    yl[ j ] = ptr_intst[ idx - 1 ] - mu;
    for( int i = 0; i < 6; ++ i ) {
      int nidx =  ptr_nidx[ 6 * ( idx - 1 ) + i ];
      if( region.count( nidx ) ) {
        yln[ i * nrow + j ] = ptr_nintst[ 6 * ( idx - 1 ) + i ] - mu;
      } else {
        yln[ i * nrow + j ] = 0;
      }
    }
  }
  
  // yl - yln * thetal
  double *vtheta = new double[ ncol ]; 
  for( int i = 0; i < 6; ++ i ) {
    vtheta[ i ] = theta[ i ];
  }
  int incx = 1;
  double alpha = - 1;
  double beta = 1;
  F77_CALL( dgemv )( "n", &nrow, &ncol, &alpha, yln, &nrow, vtheta, &incx,
            &beta, yl, &incx );

  double energy = 0;
  for( int i = 0; i < nrow; ++ i ) {
    energy += pow( yl[ i ], 2 );
  }
  energy /=  ( 2 * sigma2 );
  
  energy += log( sigma2 ) * nrow / (double)2;
  
  energy = energy + ( a + 1 ) * log( mu - mk1 ) + b / ( mu - mk1 ) - 
    log( pow( b, a ) / tgamma( a ) );
  
  energy = energy + ( alphak + 1 ) * log( sigma2 ) + betak / sigma2 -
    log( pow( betak, alphak ) / tgamma( alphak ) );
  // Rprintf( "energy = %f\n", energy );
  double pi = 3.1415926;
  double sum_theta = 0;
  if( nrow > 1) {
    for( int i = 0; i < 6; ++ i ) {
      sum_theta += pow( vtheta[ i ], 2 ); 
    }
  }
  
  energy = energy + log( 2 * pi ) * ncol / 2 + log( lambda2 ) * ncol / 2 + 
    sum_theta / ( 2 * lambda2 );

  delete [] yln;
  delete [] yl;
  delete [] vtheta;

  return energy;
}

// Single point energyY for healty voxel
double energyY( const int curr_label,
                const int curr_idx,
                double mu,
                double sigma2,
                const int *ptr_seg,
                const int *ptr_nidx,
                const double *ptr_intst,
                const double *ptr_nintst,
                const vector<double> &theta ) {
  double curr_intst = ptr_intst[ curr_idx - 1 ];
  curr_intst -= mu;
  double nbr_intst[ 6 ];
  
  for( int i = 0; i < 6; ++ i ) {
    int nidx =  ptr_nidx[ 6 * ( curr_idx - 1 ) + i ];
    if( nidx != NA_INTEGER ) {
      int nlabel = ptr_seg[ 2 * ( nidx - 1 ) ];
      if( nlabel == curr_label ) {
        nbr_intst[ i ] = ptr_nintst[ 6 * ( curr_idx - 1 ) + i ] - mu;
      } else {
        nbr_intst[ i ] = 0;
      }
    } else {
      nbr_intst[ i ] = 0;
    }
  }
  
  double energy = curr_intst;
  for( int i = 0; i < 6; ++ i ) {
    energy -= theta[ i ] * nbr_intst[ i ];
  }
  energy = pow( energy, 2 ) / ( 2 * sigma2 ) + log( sigma2 ) / 2;
  return energy;
}

// Single point energyY for tumor or outlier
double energyY( const int curr_label,
                const int curr_idx,
                double mu,
                double mk1,
                double sigma2,
                double lambda2,
                const int *ptr_seg,
                const int *ptr_nidx,
                const double *ptr_intst,
                const double *ptr_nintst,
                const vector<double> &theta,
                double alphak,
                double betak,
                double a, double b ) {
  double curr_intst = ptr_intst[ curr_idx - 1 ];
  curr_intst -= mu;
  double nbr_intst[ 6 ];
  int count = 0;
  
  for( int i = 0; i < 6; ++ i ) {
    int nidx =  ptr_nidx[ 6 * ( curr_idx - 1 ) + i ];
    if( nidx != NA_INTEGER ) {
      int nlabel = ptr_seg[ 2 * ( nidx - 1 ) ];
      if( nlabel == curr_label ) {
        nbr_intst[ i ] = ptr_nintst[ 6 * ( curr_idx - 1 ) + i ] - mu;
        ++ count;
      } else {
        nbr_intst[ i ] = 0;
      }
    } else {
      nbr_intst[ i ] = 0;
    }
  }
  
  double energy = curr_intst;
  for( int i = 0; i < 6; ++ i ) {
    energy -= theta[ i ] * nbr_intst[ i ];
  }
  energy = pow( energy, 2 ) / ( 2 * sigma2 ) + log( sigma2 ) / 2;
  
  if( count == 0 ) {
    energy = energy + ( a + 1 ) * log( mu - mk1 ) + b / ( mu - mk1 ) - 
      log( pow( b, a ) / tgamma( a ) );
    
    energy = energy + ( alphak + 1 ) * log( sigma2 ) + betak / sigma2 -
      log( pow( betak, alphak ) / tgamma( alphak ) );
    int ncol = 6;
    double pi = 3.1415926;
    double sum_theta = 0;
    for( int i = 0; i < 6; ++ i ) {
      sum_theta += pow( theta[ i ], 2 ); 
    }
    
    energy = energy + log( 2 * pi ) * ncol / 2 + log( lambda2 ) * ncol / 2 + 
      sum_theta / ( 2 * lambda2 );
  }
  
  
  return energy;
}

