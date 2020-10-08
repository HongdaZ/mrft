#include <cassert>
#include "quickhull.hpp"
#include "point.hpp"
#include <vector>

using std::vector;

// Find convex hull
vector<int> cHull( vector<int> &a ) {
  using P = point<int>;
  using S = std::vector<P>;
  using I = typename S::iterator;
  
  int len = a.size() / 2;
  S bag( len );
  
  for( int i = 0; i < len; ++ i ) {
    bag[ i ].x = a[ 2 * i ];
    bag[ i ].y = a[ 2 * i + 1 ];
  }
  I rest = quickhull::solve( bag.begin(), bag.end() );
  auto h = rest - bag.begin();
  assert( h == 2 );
  
  vector<int> res( 2 * h );
  for( int i = 0; i < h; ++ i ) {
    res[ 2 * i ] = bag[ i ].x;
    res[ 2 * i + 1 ] = bag[ i ].y;
  }
  return res;
}