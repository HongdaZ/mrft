#include <cmath>

#include "spread.h"
#include "radius.h"

using std::sqrt;

// Measure the spread of the voxels in region
double spread( const vector<int> &region, const int *ptr_aidx ) {
  
  double r = radius( region.size() );
  double max_dist = 0, dist = 0;
  const int len = region.size();
  double cr = 0, cc = 0, cs = 0, r1, c1, s1, n = 0;
  int index1;
  for( int i = 0; i < len; ++ i ) {
    if( region[ i ] != 0 ) {
      index1 = region[ i ];
      r1 = ptr_aidx[ 3 * ( index1 - 1 ) ];
      c1 = ptr_aidx[ 3 * ( index1 - 1 ) + 1 ];
      s1 = ptr_aidx[ 3 * ( index1 - 1 ) + 2 ];
      cr += r1;
      cc += c1;
      cs += s1;
      ++ n;
    }
  }
  cr /= n;
  cc /= n;
  cs /= n;
  for( int i = 0; i < len; ++ i ) {
    if( region[ i ] != 0 ) {
      index1 = region[ i ];
      r1 = ptr_aidx[ 3 * ( index1 - 1 ) ];
      c1 = ptr_aidx[ 3 * ( index1 - 1 ) + 1 ];
      s1 = ptr_aidx[ 3 * ( index1 - 1 ) + 2 ];
      dist =  sqrt( pow( r1 - cr, 2 ) + pow( c1 - cc, 2) +
        pow( s1 - cs, 2 ) );
      if( max_dist < dist ) {
        max_dist = dist;
      }
    }
  }
  
  return max_dist / r;
}