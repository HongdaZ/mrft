#include <R.h>
#include <Rinternals.h>
#include <cmath>   
#include <algorithm> 
#include <list>
#include <vector>
#include <map>

#include "updateMu.h"
#include "root.h"

using std::list;
using std::vector;
using std::map;
using std::abs;
using std::max;

// region starts from 1
// update mu for healthy cells
double updateMu( const int n,  
                 const double sigma2,
                 const double m,
                 const double nu2,
                 const double sum_theta,
                 const double sum_y ) {
  double mu = ( 1 / sigma2 * pow( 1 - sum_theta, 2 ) * sum_y + m / nu2 ) /
    ( 1 / sigma2 * n * pow( 1 - sum_theta , 2 ) + 1 / nu2 );
  return mu;
  
}

// update mu for tumor cells
double updateMu( const int n,  
                 const double sigma2,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double sum_theta,
                 const double sum_y ) {
  double y_ = sum_y / n;
  double d_ = 1 - sum_theta;
  double p;
  if( d_ == 0 ) {
    p = m;
  } else {
    double a_ = mk_1 - y_; 
    double b_ = ( a + 1 ) / ( n / sigma2 * pow( d_ , 2 ) );
    double c_ = - b / ( n / sigma2 * pow( d_ , 2 ) );
    vector<double> x( 3, 0 );
    int res = root( x, a_, b_, c_ );
    if( res == 1 ) {
      p = x[ 0 ] + mk_1;
    } else if( res == 2 ) {
      if( ( x[ 0 ] + mk_1 ) > max( y_, mk_1 ) && ( x[ 0 ] + mk_1 ) < m ) {
        p = x[ 0 ] + mk_1;
      } else {
        p = x[ 1 ] + mk_1;
      }
    } else { 
      if( ( x[ 0 ] + mk_1 ) > max( y_, mk_1 ) && ( x[ 0 ] + mk_1 ) < m ) {
        p = x[ 0 ] + mk_1;
      } else if( ( x[ 1 ] + mk_1 ) > max( y_, mk_1 ) 
                   && ( x[ 1 ] + mk_1 ) < m ) {
        p = x[ 1 ] + mk_1;
      } else {
        p = x[ 2 ] + mk_1;
      }
    }
  }
  return p;
}

// update mu for outliers
double updateMu( const int n,  
                 const double sigma2,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double sum_y ) {
  double sum_theta = 0;
  double y_ = sum_y / n;
  double d_ = 1 - sum_theta;
  double p;
  if( d_ == 0 ) {
    p = m;
  } else {
    double a_ = mk_1 - y_; 
    double b_ = ( a + 1 ) / ( n / sigma2 * pow( d_ , 2 ) );
    double c_ = - b / ( n / sigma2 * pow( d_ , 2 ) );
    vector<double> x( 3, 0 );
    int res = root( x, a_, b_, c_ );
    if( res == 1 ) {
      p = x[ 0 ] + mk_1;
    } else if( res == 2 ) {
      if( ( x[ 0 ] + mk_1 ) > max( y_, mk_1 ) && ( x[ 0 ] + mk_1 ) < m ) {
        p = x[ 0 ] + mk_1;
      } else {
        p = x[ 1 ] + mk_1;
      }
    } else { 
      if( ( x[ 0 ] + mk_1 ) > max( y_, mk_1 ) && ( x[ 0 ] + mk_1 ) < m ) {
        p = x[ 0 ] + mk_1;
      } else if( ( x[ 1 ] + mk_1 ) > max( y_, mk_1 ) 
                   && ( x[ 1 ] + mk_1 ) < m ) {
        p = x[ 1 ] + mk_1;
      } else {
        p = x[ 2 ] + mk_1;
      }
    }
  }
  return p;
}