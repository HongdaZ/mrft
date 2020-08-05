#include "getRegion.h"
#include "clearVector.h"

// get the region with label = label
void getRegion( vector<int> &region, const int &label, const int *ptr_seg,
                const int &len ){
  clearVector( region );
  for( int i = 1; i <= len; ++ i ) {
    if( ptr_seg[ 2 * ( i - 1 ) ] == label ) {
      region.push_back( i );
    }
  }
  return;
}

void getRegion( vector<int> &region, const list<int> &t_region ){
  clearVector( region );
  list<int>::const_iterator it = t_region.begin();
  ++ it;
  for( ; it != t_region.end(); ++ it ) {
    region.push_back( *it );
  }
  return;
}