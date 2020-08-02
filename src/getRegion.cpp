#include "getRegion.h"
#include "clearVector.h"

// get the region with label = label
void getRegion( vector<int> &region, const int &label,
                const int *ptr_seg, const int &len, const int row = 0 ){
  clearVector( region );
  for( int i = 1; i <= len; ++ i ) {
    if( ptr_seg[ 2 * ( i - 1 ) + row ] == label ) {
      region.push_back( i );
    }
  }
  return;
}