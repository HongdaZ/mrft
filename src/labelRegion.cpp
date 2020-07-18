#include "labelRegion.h"


// the voxel in region is labelled as 3 in ptr_seg[ 2, ]
void labelRegion( const list<int> &region, int *ptr_seg ) {
  
  for( list<int>::const_iterator it = region.begin();
       it != region.end(); ++ it ) {
    ptr_seg[ 2 * *it - 1 ] = 3;
  }
  return;
}
// change ptr_seg[ 2, region ] back
void recoverLabel( const list<int> &region, int *ptr_seg ) {
  
  for( list<int>::const_iterator it = region.begin();
       it != region.end(); ++ it ) {
    ptr_seg[ 2 * *it - 1 ] = 0;
  }
  return;
}