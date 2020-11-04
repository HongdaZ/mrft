#include <cmath>
#include "spread.h"

using std::cbrt;
using std::sqrt;

// Measure the spread of the voxels in region
double spread( const vector<int> &region, const int *ptr_aidx ) {
  const double PI = 3.141592653589793;
  double diameter = 2 * cbrt( ( double )region.size() * 3 / 4 / PI  );
  double max_dist = 0, dist = 0;
  const int len = region.size();
  int r1, c1, s1, r2, c2, s2;
  int index1, index2;
  for( int i = 0; i < len - 1; ++ i ) {
    if( region[ i ] != 0 ) {
      index1 = region[ i ];
      r1 = ptr_aidx[ 3 * ( index1 - 1 ) ];
      c1 = ptr_aidx[ 3 * ( index1 - 1 ) + 1 ];
      s1 = ptr_aidx[ 3 * ( index1 - 1 ) + 2 ];
      for( int j = i + 1; j < len; ++ j ) {
        if( region[ j ] != 0 ) {
          index2 = region[ j ];
          r2 = ptr_aidx[ 3 * ( index2 - 1 ) ];
          c2 = ptr_aidx[ 3 * ( index2 - 1 ) + 1 ];
          s2 = ptr_aidx[ 3 * ( index2 - 1 ) + 2 ];
          dist = sqrt( pow( r1 - r2, 2 ) + pow( c1 - c2, 2) +
            pow( s1 - s2, 2 ) );
          if( max_dist < dist ) {
            max_dist = dist;
          }
        }
      }
    }
  }
  return max_dist / diameter;
}