#include <cmath>

#include "roundness.h"
#include "radius.h"

using std::sqrt;

// Measure the roundness of the voxels in region
double roundness( const vector<int> &region, const int *ptr_aidx ) {
  
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
      dist =  sqrt( pow( r1 - cr + 0.71, 2 ) + pow( c1 - cc + 0.71, 2 ) +
        pow( s1 - cs + 0.71, 2 ) );
      if( max_dist < dist ) {
        max_dist = dist;
      }
    }
  }
  
  max_dist /= r;
  return max_dist * max_dist;
}
// 2D version of the function above
double roundness( const vector<int> &region, const int &plane,
               const int *ptr_aidx ) {
  vector<int> plane_idx{ 1, 2, 0, 2, 0, 1 };
  vector<int> curr_plane{ plane_idx[ 2 * plane ], 
                          plane_idx[ 2 * plane + 1 ] };
  double r = radius2D( region.size() );
  double max_dist = 0, dist = 0;
  const int len = region.size();
  double n = 0;
  vector<double> center( 2, 0 );
  vector<double> xyz( 2, 0 );
  int index1;
  for( int i = 0; i < len; ++ i ) {
    if( region[ i ] != 0 ) {
      index1 = region[ i ];
      for( int j = 0; j < 2; ++ j ) {
        xyz[ j ] = ptr_aidx[ 3 * ( index1 - 1 ) + curr_plane[ j ] ];
        center[ j ] += xyz[ j ];
      }
      ++ n;
    }
  }
  center[ 0 ] /= n;
  center[ 1 ] /= n;
  
  for( int i = 0; i < len; ++ i ) {
    if( region[ i ] != 0 ) {
      index1 = region[ i ];
      for( int j = 0; j < 2; ++ j ) {
        xyz[ j ] = ptr_aidx[ 3 * ( index1 - 1 ) + curr_plane[ j ] ];
      }
      dist =  sqrt( pow( xyz[ 0 ] - center[ 0 ] + 0.71, 2 ) + 
        pow( xyz[ 1 ] - center[ 1 ] + 0.71, 2 ) );
      if( max_dist < dist ) {
        max_dist = dist;
      }
    }
  }
  max_dist /= r;
  return max_dist * max_dist;
}