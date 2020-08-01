#include "assignRegion.h"

// store region in regions
void assignRegion( vector<int> &regions, vector<int> &region, 
                   const int &row, int label ) {
  for( int i = 0; i < region.size(); ++ i ) {
    regions[ 2 * ( region[ i ] - 1 ) + row ] = label;
  }
  return;
}