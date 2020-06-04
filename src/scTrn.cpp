#include <R.h>
#include <Rinternals.h>

#include "scTrn.h"

// start = 1 to length array
int scTrn( list< map<int, int>> &regions, 
           int *ptr_label, const int *ptr_nidx, int start ) {
  int current = ptr_label[ 2 * ( start - 1 ) ];
  if( current <= - 1 && current >= -3 ) {
    return 0;
  } else if( current <= - 4 ) {
    // split
  }
  
}