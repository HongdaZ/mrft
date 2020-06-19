#include <R.h>
#include <Rinternals.h>

#include <cmath>   

#include "updateParm.h"

using std::abs;

// update parameters for healthy cells
void updateParm( double &mu, vector<double> &theta, double &sigma2, 
                 set<int> &region,
                 double m,
                 double nu2,
                 const double *ptr_intst,
                 int curr_label,
                 double lambda2,
                 const int *ptr_seg,
                 const int *ptr_nidx,
                 const double *ptr_nintst,
                 double alphal,
                 double betal,
                 int maxit ) {

  for( int i = 0; i < 6; ++ i ) {
    theta[ i ] = 0;
  }
  int i = 0;
  double tol = 1;
  double tmp = 0;
  mu = - 1;
  while( i < maxit && tol > .0001 ) {
    tmp = updateMu( region, sigma2, m, nu2, theta, ptr_intst );
    // Rprintf( "mu = %f \n", tmp );
    tol = abs( mu - tmp );
    mu = tmp;
    updateTS( region, curr_label, mu, sigma2, lambda2, ptr_seg, ptr_nidx,
              ptr_intst, ptr_nintst, theta, alphal, betal );
    // for( int j = 0; j < 6; ++ j ) {
    //   Rprintf( "%f\t", theta[ j ] );
    // }
    // Rprintf( "\n" );
    // Rprintf( "sigma2 = %f\n", sigma2 );
    ++ i;
  }
  
  
}