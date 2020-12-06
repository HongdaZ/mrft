#include "pTNbr.h"

#include "Rinternals.h"

// Find the proportion of new voxels having tumor neighbor
double pTNbr( const vector<int> &region, const int *ptr_seg1, 
           const int &label1, const int *ptr_nidx ) {
  int total = 0, n = 0, idx = 0, n_idx = 0;
  for( vector<int>::const_iterator it = region.begin();
       it != region.end(); ++ it ) {
    idx = *it;
    if( idx != 0 ) {
      ++ total;
      for( int i = 0; i < 6; ++ i ) {
        n_idx = ptr_nidx[ 6 * ( idx - 1 ) + i ];
        if( n_idx !=  NA_INTEGER ) {
          if( ptr_seg1[ 2 * ( n_idx - 1 ) ] == label1 ) {
            ++ n;
            break;
          }
        }
      } 
    }
  }
  return (double)n / total;
}