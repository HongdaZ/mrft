#include "skip.h"

// whether skip current voxel
bool skip( const int idx, const int *ptr_seg, const int *ptr_nidx ) {
  int curr_label = ptr_seg[ 2 * ( idx - 1 ) ];
  list<int> nbr_label;
  list<int> tumor_nbr;
  nbrLabel( nbr_label, tumor_nbr, idx, ptr_seg, ptr_nidx );
  bool res = true;
  if( nbr_label.size() > 0 ) {
    for( list<int>::iterator it = nbr_label.begin(); 
         it != nbr_label.end(); ++ it ) {
      if( *it != 0 ) {
        res = false;
      }
    }
  }
  
  res = res && ( curr_label == 0 );
  return res;
}