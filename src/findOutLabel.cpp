#include "findOutLabel.h"

// find an outlier label for current voxel
int findOutLabel( const int idx, const int *ptr_seg, 
                  const list<int> &outl_labels ) {
  int out_label;
  // current label is not outlier
  if( ptr_seg[ 2 * ( idx - 1 ) ] < 1 ) {
    if( outl_labels.empty() ) {
      out_label = 1;
    } else if( *( -- outl_labels.end() ) == outl_labels.size() ){
      out_label = outl_labels.size() + 1;
    } else {
      int i = 0;
      list<int>::const_iterator it = outl_labels.begin();
      for( ; it != outl_labels.end(); ++ it, ++ i ) {
        if( ( i + 1 ) < *it ) {
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