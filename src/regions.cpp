#include "regions.h"
#include "cnctRegion.h"

// Find all the connected regions in ptr_seg with label == label
list<vector<int>> regions( int *ptr_seg, const int &len, 
                           vector<int> &region,
                           const int &label,
                           const int *ptr_nidx,
                           const int *ptr_aidx ) {
  list<vector<int>> res;
  int n;
  int idx;
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_seg, ptr_seg, 
                    label, region ) ) {
      n = region.size();
      vector<int> r( n * 4 );
      for( int j = 0; j < n; ++ j ) {
        idx = region[ j ];
        r[ 4 * j ] = idx;
        for( int k = 0; k < 3; ++ k ) {
          r[ 4 * j + k + 1 ] = ptr_aidx[ 3 * ( idx - 1 ) + k ];
        }
      }
      res.push_back( r );
    }
  }
  return res;
}