#include "Rinternals.h"

#include "peel.h"

int peel( int *ptr_whole, const int &len,
           const int *ptr_nidx ) {
  int last = -1, n_tumor = 1, nidx = 0;
  bool rm = false;
  // Remove all tumor voxels
  while( n_tumor > 0 ) {
    n_tumor = 0;
    for( int i = 0; i < len; ++ i ) {
      if( ptr_whole[ 2 * i ] == 1 ) {
        rm = false;
        for( int j = 0; j < 6; ++ j ) {
          nidx = ptr_nidx[ 6 * i + j ];
          if( nidx != NA_INTEGER ) {
            if( ptr_whole[ 2 * ( nidx - 1 ) ] == 0 ) {
              rm = true;
              ptr_whole[ 2 * i ] = 0;
              last = i;
              break;
            }
          } else {
            rm = true;
            ptr_whole[ 2 * i ] = 0;
            last = i;
            break;
          }
        }
        if( ! rm ) {
          ++ n_tumor;
        }
      }
    }
  }
  return last;
}