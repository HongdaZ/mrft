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
        if( ptr_seg1[ 2 * ( n_idx - 1 ) ] == 0 &&
            ptr_seg2[ 2 * ( n_idx - 1 ) ] == label ) {
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
void excldRegion( const vector<int> &region,
                  int *ptr_seg1,
                  const int *ptr_seg2, 
                  const int &label, const int &size ) {
  int index, n = 0;
  for( vector<int>::const_iterator it = region.begin(); 
       it != region.end(); ++ it ) {
    index = *it;
    if( ptr_seg2[ 2 * ( index - 1 ) ] == label ) {
      ++ n;
      if( n > size ) {
        return;
      }
    }
  }
  for( vector<int>::const_iterator it = region.begin(); 
       it != region.end(); ++ it ) {
    index = *it;
    ptr_seg1[ 2 * ( index - 1 ) ] = 0;
  }
}

void excldRegion( const vector<int> &region,
                  int *ptr_seg, const int &size ) {
  if( region.size() > size ) {
    return;
  } 
  int index;
  for( vector<int>::const_iterator it = region.begin(); 
       it != region.end(); ++ it ) {
    index = *it;
    ptr_seg[ 2 * ( index - 1 ) ] = 0;
  }
}