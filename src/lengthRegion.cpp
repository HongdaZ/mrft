#include "lengthRegion.h"
// get the length of existing region 
int lengthRegion( const int *ptr_seg, const int &len,
                  const int curr_label ) {
  int len_region = 0;
  for( int k = 0; k < len; ++ k ) {
    if( ptr_seg[ 2 * k ] == curr_label ) {
      ++ len_region; 
    }
  }
  return len_region;
}