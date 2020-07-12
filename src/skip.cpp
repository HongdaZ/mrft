#include "skip.h"

// whether skip current voxel
bool skip( const int idx, const int *ptr_seg, const int *ptr_nidx ) {
  int curr_label = ptr_seg[ 2 * ( idx - 1 ) ];
  vector<int> nbr_label = nbrLabel( idx, ptr_seg, ptr_nidx );
  bool res = true;
  for( int i = 0; i < 6; ++ i ) {
    if( nbr_label[ i ] != NA_INTEGER && nbr_label[ i ] != 0 ) {
      res = false;
    }
  }
  res = res && ( curr_label == 0 );
  return res;
}