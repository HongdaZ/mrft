#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>

#include <set>
#include <vector>

#include "updateTheta.h"

using std::set;
using std::vector;



// region starts from 1
// updateTheta for health regions
vector<double> updateTheta( set<int> &region,
                            double mu,
                            double sigma2,
                            double lambda2,
                            const double *ptr_seg,
                            const double *ptr_nidx,
                            const double *ptr_intst,
                            const double *ptr_nintst ) {
  int nrow = region.size();
  int ncol = 6;
  double *yln = new double[ nrow * ncol ];
  double *yl = new double[ nrow ];
  int j = 0;
  set<int>::iterator it = region.begin();
  int idx = *it;
  int label = ptr_seg[ 2 * ( idx - 1 ) ];
  // Initialize yln and yl;
  for( ; it!= region.end(); ++ it, ++ j ) {
    idx = *it;
    yl[ j ] = ptr_intst[ idx - 1 ] - mu;
    for( int i = 0; i < 6; ++ i ) {
      double tmp = ptr_nintst[ 6 * ( idx - 1 ) + i ];
      if( ISNAN( tmp ) ) {
        tmp = 0;
      } else {
        int nidx =  ptr_nidx[ 6 * ( idx - 1 ) + i ];
        if( nidx != NA_INTEGER ) {
          int nlabel = ptr_seg[ 2 * ( nidx - 1 ) ];
          if( nlabel == label ) {
            yln[ i * nrow + j ] = tmp;
          } else {
            yln[ i * nrow + j ] = 0;
          }
        }
      }
      
    }
    
  }
  
  delete [] yln;
  delete [] yl;
  
}
