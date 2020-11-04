#include <R.h>
#include <Rinternals.h>
#include "onRegion.h"
#include "cnctRegion.h"
#include "excldRegion.h"
#include "sliceCount.h"
#include "zeroVector.h"
#include "pad2zero.h"
#include "spread.h"

// Extend ptr_seg1 to adjacent regions of ptr_seg2
void onRegion( int *ptr_on, const int &len, const double &prop,
               int *ptr_seg1, const int &label1, 
               int *ptr_seg2, const int &label2,
               vector<int> &region, 
               const int *ptr_nidx, const int *ptr_aidx,
               const int &nr, const int &nc, const int &ns ) {
  vector<int> bound{ 0, nr - 1, nr + nc - 1, nr + nc + ns - 1 };
  vector<int> seg1_count( nr + nc + ns, 0 );
  vector<int> seg2_count( nr + nc + ns, 0 );
  sliceCount( seg1_count, ptr_seg1, len, label1, ptr_nidx, ptr_aidx,
              nr, nc, ns );
  int n_vio = 0, idx = 1;
  double spread_idx = 0;
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_seg2, ptr_seg2,
                    label2, region ) ) {
      if( ! excldRegion( region, ptr_nidx, ptr_seg2,
                         ptr_seg1, label1 ) ) {
        n_vio = 0;
        zeroVector( seg2_count );
        sliceCount( seg2_count, region, ptr_nidx, ptr_aidx, 
                    nr, nc, ns );
        for( int j = 0; j < 3; ++ j ) {
          for( int k = bound[ j ]; k < bound[ j + 1 ]; ++ k ) {
            if( seg2_count[ k ] > prop * seg1_count[ k ] ) {
              ++ n_vio;
              break;
            }
          }
        }
        
        if( n_vio < 2 ) {
          spread_idx = spread( region, ptr_aidx );
          if( spread_idx < 3 ) {
            for( vector<int>::const_iterator it = region.begin();
                 it != region.end(); ++ it ) {
              idx = *it;
              if( idx != 0 ) {
                ptr_on[ 2 * ( idx - 1 ) ] = 1;
              }
            }
          }
        }
      }
    }
  }
  pad2zero( ptr_seg2, len );
}