#include "pad2zero.h"

// Recover the padding to zero's
void pad2zero( int *ptr_seg, const int &len ) {
  for( int i = 0; i < len; ++ i ) {
    ptr_seg[ 2 * i + 1 ] = 0;
  }
}