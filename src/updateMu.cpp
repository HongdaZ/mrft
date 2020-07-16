#include <R.h>
#include <Rinternals.h>
#include <cmath>   

#include <list>
#include <vector>
#include <map>

#include "updateMu.h"

using std::list;
using std::vector;
using std::map;
using std::abs;

// region starts from 1
// update mu for healthy cells
double updateMu( const list<int> &region,  
                 const double sigma2,
                 const double m,
                 const double nu2,
                 const vector<double> theta,
                 const double *ptr_intst ) {
  int n = region.size();
  double sum_theta = 0;
  for( vector<double>::const_iterator it = theta.begin(); it != theta.end(); 
    ++ it ) {
    sum_theta += *it;
  }
  double sum_y = 0;
  for( list<int>::const_iterator it = region.begin(); it != region.end();
  ++ it ) {
    sum_y += ptr_intst[ *it - 1 ]; 
  }
  double mu = ( 1 / sigma2 * pow( 1 - sum_theta, 2 ) * sum_y + m / nu2 ) /
    ( 1 / sigma2 * n * pow( 1 - sum_theta , 2 ) + 1 / nu2 );
  return mu;
  
}
double derivative( double mu, double sigma2, double sum_theta, 
             double mean_y, int n, double a, double b, double m, double mk ) {
  return n / sigma2 * pow( 1 - sum_theta, 2 ) * ( mean_y - mu ) -
    ( a + 1 ) / ( mu - m ) + b / pow( mu - m, 2 );
}
// update mu for tumor cells
double updateMu( const list<int> &region,  
                 const double sigma2,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const vector<double> theta,
                 const double *ptr_intst ) {
  int n = region.size();
  double sum_theta = 0;
  for( vector<double>::const_iterator it = theta.begin(); it != theta.end(); 
  ++ it ) {
    sum_theta += *it;
  }
  double sum_y = 0;
  for( list<int>::const_iterator it = region.begin(); it != region.end();
  ++ it ) {
    sum_y += ptr_intst[ *it - 1 ]; 
  }
  int max_itr = 100;
  double l, r;
  double mean_y = sum_y / n;
  if( mean_y > mk_1 ) {
    l = mean_y;
  } else {
    l = mk_1 + (  m- mk_1 ) / 1000;
  }
  r = m;
  int i = 0;
  double dif =  r - l;
  double p = ( l + r ) / 2;
  double tol = (  m- mk_1 ) / 1000 ;
  
  // Rprintf( "mean_y = %f \n", mean_y );
  double fr = derivative( r, sigma2, sum_theta, mean_y, n, a, b, mk_1,  m);
  double fl = derivative( l, sigma2, sum_theta, mean_y, n, a, b, mk_1,  m);
  if( fl < 0 && fr < 0 ) {
    return l;
  }
  double fp;
  while( abs( dif ) > tol && i < max_itr &&
         derivative( p, sigma2, sum_theta, mean_y, n, a, b, mk_1, m ) != 0 ) {
    fp = derivative( p, sigma2, sum_theta, mean_y, n, a, b, mk_1, m );
    fl = derivative( l, sigma2, sum_theta, mean_y, n, a, b, mk_1, m );
    if( fp * fl < 0 ) {
      r = p;
    } else {
      l = p;
    }
    dif = r - l;
    p = ( r + l ) / 2;
    ++ i;
    // Rprintf( "fp= %f\n", fp );
  }
  return p;
}

// update mu for outliers
double updateMu( const list<int> &region,  
                 const double sigma2,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double *ptr_intst ) {
  int n = region.size();
  double sum_theta = 0;
  double sum_y = 0;
  for( list<int>::const_iterator it = region.begin(); it != region.end();
  ++ it ) {
    sum_y += ptr_intst[ *it - 1 ]; 
  }
  int max_itr = 100;
  double l, r;
  double mean_y = sum_y / n;
  if( mean_y > mk_1 ) {
    l = mean_y;
  } else {
    l = mk_1 + ( m - mk_1 ) / 1000;
  }
  r = m;
  int i = 0;
  double dif =  r - l;
  double p = ( l + r ) / 2;
  double tol = ( m - mk_1 ) / 1000 ;
  
  // Rprintf( "mean_y = %f \n", mean_y );
  double fr = derivative( r, sigma2, sum_theta, mean_y, n, a, b, mk_1, m );
  double fl = derivative( l, sigma2, sum_theta, mean_y, n, a, b, mk_1, m );
  if( fl < 0 && fr < 0 ) {
    return l;
  }
  double fp;
  while( abs( dif ) > tol && i < max_itr &&
         derivative( p, sigma2, sum_theta, mean_y, n, a, b, mk_1, m ) != 0 ) {
    fp = derivative( p, sigma2, sum_theta, mean_y, n, a, b, mk_1, m );
    fl = derivative( l, sigma2, sum_theta, mean_y, n, a, b, mk_1, m );
    if( fp * fl < 0 ) {
      r = p;
    } else {
      l = p;
    }
    dif = r - l;
    p = ( r + l ) / 2;
    ++ i;
    // Rprintf( "fp= %f\n", fp );
  }
  return p;
}