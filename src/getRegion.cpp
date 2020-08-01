#include "getRegion.h"

// get the region with label = label
void getRegion( vector<int> &region, const int &len_region, const int &label,
                const int *ptr_seg, const int &len ){
  int count = 0;
  for( int i = 1; i <= len; ++ i ) {
    if( ptr_seg[ 2 * ( i - 1 ) ] == label ) {
      region[ count ] = i;
      ++ count;
      if( count == len_region ) {
        break;
      }
    }
  }
  return;
}