#include "findOutLabel.h"

// find an outlier label for current voxel
int findOutLabel( const int idx, const int *ptr_seg, 
                  vector<int> &outl_labels ) {
  int out_label;
  // current label is not outlier
  if( ptr_seg[ 2 * ( idx - 1 ) ] < 1 ) {
    for( int i = 0; i < outl_labels.size(); ++ i ) {
      if( outl_labels[ i ] == 0 ) {
        out_label = i + 1;
      }
    }
  } else {
    out_label = ptr_seg[ 2 * ( idx - 1 ) ];
  }
  return out_label;
}