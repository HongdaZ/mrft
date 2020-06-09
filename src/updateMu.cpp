#include <R.h>
#include <Rinternals.h>

#include <set>
#include <vector>

#include "updateMu.h"

using std::set;
using std::vector;

// update mu for healthy cells
double updateMu( set<int> &region,  
           double sigma2,
           double m,
           double nu2,
           vector<double> theta,
           double *ptr_intst ) {
  int n = region.size();
  double sum_theta = 0;
  for( vector<double>::iterator it = theta.begin(); it != theta.end(); 
    ++ it ) {
    sum_theta += *it;
  }
  double sum_y = 0;
  for( set<int>::iterator it = region.begin(); it != region.end();
  ++ it ) {
    sum_y += ptr_intst[ *it - 1 ]; 
  }
  double mu = ( 1 / sigma2 * pow( 1 - sum_theta, 2 ) * sum_y + m / nu2 ) /
    ( 1 / sigma2 * n * pow( 1 - sum_theta , 2 ) + 1 / nu2 );
  return mu;
  
}