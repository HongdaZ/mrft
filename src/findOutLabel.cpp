#include "findOutLabel.h"

// find an outlier label for current voxel
int findOutLabel( const int idx, const int *ptr_seg, 
                 const set<int> &outl_labels ) {
  int out_label;
  // current label is not outlier
  if( ptr_seg[ 2 * ( idx - 1 ) ] < 1 ) {
    if( outl_labels.empty() ) {
      out_label = 1;
    } else if( *( -- outl_labels.end() ) == outl_labels.size() ){
      out_label = outl_labels.size() + 1;
    } else {
      int i = 0;
      set<int>::iterator it_set = outl_labels.begin();
      for( ; it_set != outl_labels.end(); ++ it_set, ++ i ) {
        if( ( i + 1 ) < *it_set ) {
          out_label = i + 1;
          break;
        }
      }
    }
  } else {
    out_label = ptr_seg[ 2 * ( idx - 1 ) ];
  }
  return out_label;
}