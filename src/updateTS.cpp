#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include <algorithm>

#include <list>
#include <vector>

#include "updateTS.h"

using std::list;
using std::vector;
using std::find;


// region starts from 1
// updateTheta and sigma2 for health and tumorous regions
void updateTS( const list<int> &region, int curr_label,
               const double mu,
               double &sigma2,
               const double lambda2,
               const int *ptr_seg,
               const int *ptr_nidx,
               const double *ptr_intst,
               const double *ptr_nintst,
               vector<double> &theta,
               const double alphal,
               const double betal ){
  int nrow = region.size();
  int ncol = 6;
  double *yln = new double[ nrow * ncol ];
  double *yl = new double[ nrow ];
  
  list<int>::const_iterator it = region.begin();
  int idx = *it;
  int label = ptr_seg[ 2 * ( idx - 1 ) ];
  // Initialize yln and yl;
  for( int j = 0; it!= region.end(); ++ it, ++ j ) {
    idx = *it;
    yl[ j ] = ptr_intst[ idx - 1 ] - mu;
    for( int i = 0; i < 6; ++ i ) {
      // healthy
      int nidx =  ptr_nidx[ 6 * ( idx - 1 ) + i ];
      if( curr_label < 0 && curr_label > - 4 ) {
        if( nidx != NA_INTEGER ) {
          int nlabel = ptr_seg[ 2 * ( nidx - 1 ) ];
          if( nlabel == label ) {
            yln[ i * nrow + j ] = ptr_nintst[ 6 * ( idx - 1 ) + i ] - mu;
          } else {
            yln[ i * nrow + j ] = 0;
          }
        } else {
          yln[ i * nrow + j ] = 0;
        }
        // tumor
      } else {
        list<int>::const_iterator it = find( region.begin(), 
                                             region.end(), nidx );
        if( it != region.end() ) {
          yln[ i * nrow + j ] = ptr_nintst[ 6 * ( idx - 1 ) + i ] - mu;
        } else {
          yln[ i * nrow + j ] = 0;
        }
      }
      
    }
  }
  // for( int i = 0; i < nrow; ++ i ) {
  //   Rprintf( "%f\n", yl[ i ] );
  // }
  // for( int i = 0; i < nrow; ++ i ) {
  //   Rprintf( "\n" );
  //   for( int j = 0; j < 6; ++ j ) {
  //     Rprintf( "%f,\t", yln[ nrow * j + i ] );
  //   }
  // }
  
  // Identity matrix
  double *I = new double[ ncol * ncol ];
  for( int i = 0; i < ncol; ++ i ) {
    for( int j = 0; j < ncol; ++ j ) {
      I[ i * ncol + j ] = ( i == j ) ? 1 : 0;
    }
  }
  double *tmp = new double[ ncol ];
  double *vtheta = new double[ ncol ]; 
  // t( Yln ) Yln / sigma2 + I / lambda2 
  double alpha = 1 / sigma2;
  double beta = 1 / lambda2;
  // Rprintf( "nrow = %d, ncol = %d, alpha = %f, beta = %f\n", 
  //          nrow, ncol, alpha, beta );
  F77_CALL( dsyrk )( "u", "t", &ncol, &nrow, &alpha, yln, &nrow, &beta, 
            I, &ncol );
  for( int i = 0; i < ncol; ++ i ) {
    for( int j = i + 1; j < ncol; ++ j ) {
      I[ i * ncol + j ] = I[ j * ncol + i ];
    }
  }
  
  // inverse
  int *ipiv = new int[ ncol ];
  int lwork = ncol * ncol;
  double *work = new double[ lwork ];
  int info;
  F77_CALL( dgetrf )( &ncol, &ncol, I, &ncol, ipiv, &info );
  F77_CALL( dgetri )( &ncol, I, &ncol, ipiv, work, &lwork, &info );
  // // Debug for inverse 
  // for( int i = 0; i < ncol; ++ i ) {
  //   for( int j = 0; j < ncol; ++ j ) {
  //     Rprintf( "%f\t", I[ i * ncol + j ] );
  //   }
  //   Rprintf( "\n" );
  // }
  // t( yln ) yl \ sigma2
  int incx = 1;
  double *x = new double[ nrow ];
  beta = 0;
  F77_CALL( dgemv )( "t", &nrow, &ncol, &alpha, yln, &nrow, yl, &incx, &beta, 
            tmp, &incx );
  // // debug matrix vector multiply
  // for( int i = 0; i < ncol; ++ i ) {
  //   Rprintf( "%f\t", tmp[ i ] );
  // }
  // ////////////////////////////////////
  alpha = 1;
  F77_CALL( dgemv )( "t", &ncol, &ncol, &alpha, I, &ncol, tmp, &incx, &beta, 
            vtheta, &incx );
  for( int i = 0; i < ncol; ++ i ) {
    theta[ i ] = vtheta[ i ];
  }
  // // debug matrix vector multiply
  // for( int i = 0; i < ncol; ++ i ) {
  //   Rprintf( "%f\t", theta[ i ] );
  // }
  // ////////////////////////////////////
  
  // update sigma2
  // yl - yln * thetal
  alpha = - 1;
  beta = 1;
  F77_CALL( dgemv )( "n", &nrow, &ncol, &alpha, yln, &nrow, vtheta, &incx, 
            &beta, yl, &incx );
  // // debug matrix vector multiply
  // for( int i = 0; i < nrow; ++ i ) {
  //   Rprintf( "%f\t", yl[ i ] );
  // }
  // Rprintf( "\n");
  // ////////////////////////////////////
  double sum = betal;
  for( int i = 0; i < nrow; ++ i ) {
    sum += pow( yl[ i ], 2 ) / 2;
  }
  sigma2 = sum / ( nrow / (double)2 + alphal + 1 );

  delete [] yln;
  delete [] yl;
  delete [] I;
  delete [] tmp;
  delete [] vtheta;
  delete [] ipiv;
  delete [] work;
  delete [] x;
  
  
  return;
  
}
