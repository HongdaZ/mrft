#include "extRegion.h"

// Extend the region
void extRegion( const vector<int> &region, int *ptr_seg,
                const int &label ) {
  int n1 = 0;
  int n2 = 0;
  
  for( vector<int>::const_iterator it = region.begin();
       it != region.end(); ++ it ) {
    if( *it != 0 ) {
      ++ n1;
      if( ptr_seg[ 2 * ( *it - 1 ) ] > 0 ) {
        ++ n2;
      }
    }
  }
  if( n2 > n1 * .5 ) {
    for( vector<int>::const_iterator it = region.begin();
         it != region.end(); ++ it ) {
      if( *it != 0 ) {
        ptr_seg[ 2 * ( *it - 1 ) ] = label;
      }
    }
  }
}