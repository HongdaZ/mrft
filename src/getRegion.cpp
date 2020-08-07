#include "getRegion.h"
#include "clearVector.h"

// get region and region_label
void getRegion( int &region_label, vector<int> &region, 
                const vector<int> &regions_whole, 
                const vector<int> &regions_sub,
                const int &order ) {
  clearVector( region );
  if( order == 0 ) {
    region_label = regions_whole.front();
    for( int i = 1; i < ( regions_whole.size() - 1 ); ++ i ) {
      region.push_back( regions_whole[ i ] );
    }
  } else {
    int i = 0;
    int count = 0;
    for( ; i < regions_sub.size(); ++ i ) {
      if( regions_sub[ i ] != NA_INTEGER && regions_sub[ i ] < 0 ) {
        ++ count;
        if( count == order ) {
          break;
        }
      }
    }
    region_label = regions_sub[ i ];
    ++ i;
    while( regions_sub[ i ] != NA_INTEGER ) {
      region.push_back( regions_sub[ i ] );
    }
  }
  return;
}