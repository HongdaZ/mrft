#include "excldRegion.h"

#include "Rinternals.h"

void excldRegion( const vector<int> &region, const int *ptr_nidx,
                  int *ptr_seg1,
                  const int *ptr_seg2, const int &label ) {
  int n_idx, index;
  for( vector<int>::const_iterator it = region.begin(); 
       it != region.end(); ++ it ) {
    index = *it;
    for( int i = 0; i < 6; ++ i ) {
      n_idx = ptr_nidx[ 6 * ( index - 1 ) + i ];
      if( n_idx != NA_INTEGER ) {
        if( ptr_seg2[ 2 * ( n_idx - 1 ) ] == label ) {
          return;
        }
      }
    }
  }
  for( vector<int>::const_iterator it = region.begin(); 
       it != region.end(); ++ it ) {
    index = *it;
    ptr_seg1[ 2 * ( index - 1 ) ] = 0;
  }
}