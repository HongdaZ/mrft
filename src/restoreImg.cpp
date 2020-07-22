#include <R.h>
#include <Rinternals.h>

#include "restoreImg.h"

// restore 3D image from vector
void restoreImg( const int *ptr_idx, int *ptr_res_seg, int *ptr_res_image, 
                 int len ) {
  int label, vidx;
  for( int i = 0; i < 240 * 240 * 155; ++ i ) {
    ptr_res_image[ i ] = NA_INTEGER;
  }
  for( int i = 1; i <= len; ++ i ) {
    label = ptr_res_seg[ 2 * ( i - 1 ) ];
    vidx = ptr_idx[ i - 1 ];
    ptr_res_image[ vidx - 1 ] = label;
  }
  return;
}