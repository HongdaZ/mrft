#include <list>

#include "inRegion.h"
#include "regions.h"
#include "region2slice.h"
#include "enclose.h"

using std::list;

// Find the points of ptr_seg2 enclosed in ptr_seg1
void inRegion( int *ptr_enclose, const int &len,
               int *ptr_seg1, const int &label1, 
               int *ptr_seg2, const int &label2,
               vector<int> &region, 
               const int *ptr_nidx, const int *ptr_aidx,
               const int &nr, const int &nc, const int &ns ) {
  // Find connected regions in ptr_seg1
  list<vector<int>> regions1 = regions( ptr_seg1, len,
                                        region, label1, 
                                        ptr_nidx,
                                        ptr_aidx );
  vector<list<vector<int>>> slices1 =
    region2slice( regions1, nr, nc, ns );
  list<vector<int>> regions2 = regions( ptr_seg2, len,
                                        label2, ptr_aidx ); 
  vector<list<vector<int>>> slices2 = 
    region2slice( regions2, nr, nc, ns );
  enclose( ptr_enclose, len, slices1, slices2 );
}

