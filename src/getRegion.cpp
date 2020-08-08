#include "getRegion.h"
#include "clearVector.h"

// get region and region_label
void getRegion( int &region_label, vector<int> &region, 
                const vector<int> &regions_whole, 
                const vector<int> &regions_sub,
                int &start) {
  clearVector( region );
  if( start == -1 ) {
    region_label = regions_whole.front();
    for( int i = 1; i < ( regions_whole.size() - 1 ); ++ i ) {
      region.push_back( regions_whole[ i ] );
    }
    start = 0;
  } else {
    int i = start;
    region_label = regions_sub[ i ];
    while( regions_sub[ ++ i ] != NA_INTEGER ) {
      region.push_back( regions_sub[ i ] );
    }
    start = i + 1;
  }
  return;
}