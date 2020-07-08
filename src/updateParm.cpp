#include <R.h>
#include <Rinternals.h>

#include <cmath>   

#include "updateParm.h"

using std::abs;

// update parameters for healthy cells
void updateParm( double &mu, vector<double> &theta, double &sigma2, 
                 const set<int> &region,
                 const double m,
                 const double nu2,
                 const double *ptr_intst,
                 int curr_label,
                 const double lambda2,
                 const int *ptr_seg,
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
  while( i < 2 || ( i < maxit && tol > .0001 ) ) {
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

// update parameters for tumor cells
void updateParm( double &mu, vector<double> &theta, double &sigma2, 
                 const map<int, int> &region,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double *ptr_intst,
                 const int curr_label,
                 const double lambda2,
                 const int *ptr_seg,
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
  set<int> set_region;
  for( map<int, int>::const_iterator it = region.begin(); 
       it != region.end(); ++ it ) {
    set_region.insert( it->first );
  }
  while(  i < 2 || ( i < maxit && tol > .0001 ) ) {
    tmp = updateMu( region, sigma2, m, mk_1, a, b, theta, ptr_intst );
    // Rprintf( "mu = %f \n", tmp );
    tol = abs( mu - tmp );
    mu = tmp;
    updateTS( set_region, curr_label, mu, sigma2, lambda2, ptr_seg, ptr_nidx,
              ptr_intst, ptr_nintst, theta, alphal, betal );
    // for( int j = 0; j < 6; ++ j ) {
    //   Rprintf( "%f\t", theta[ j ] );
    // }
    // Rprintf( "\n" );
    // Rprintf( "sigma2 = %f\n", sigma2 );
    ++ i;
  }
}

// update parameters for outliers
void updateParm( double &mu, double &sigma2, 
                 const map<int, int> &region,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double *ptr_intst,
                 int curr_label,
                 const double lambda2,
                 const int *ptr_seg,
                 const int *ptr_nidx,
                 const double alphal,
                 const double betal,
                 int maxit ) {
  int i = 0;
  double tol = 1;
  double tmp = 0;
  mu = - 1;
  map<int, int>::const_iterator it = region.begin();
  int idx = it->first;
  while(  i < 2 || ( i < maxit && tol > .0001 ) ) {
    tmp = updateMu( region, sigma2, m, mk_1, a, b, ptr_intst );
    // Rprintf( "mu = %f \n", tmp );
    tol = abs( mu - tmp );
    mu = tmp;
    // beta is 2.5 times as large as beta for tumor region;
    updateSigma( idx, mu, ptr_intst, alphal, 2.5 * betal, sigma2 );
    // Rprintf( "\n" );
    // Rprintf( "sigma2 = %f\n", sigma2 );
    ++ i;
  }
}