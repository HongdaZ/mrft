#include "Rinternals.h"

#include <list>

#include "peel.h"
#include "clearVector.h"
#include "zeroVector.h"

using std::list;

int peel( int *ptr_one, const int &len,
           const int *ptr_nidx ) {
  list<int> outer_old;
  list<int> outer_new;
  int last = 0, idx = 0, nidx = 0;
  bool rm = false;
  // Find the initial outer tumor voxels
  for( int i = 0; i < len; ++ i ) {
    if( ptr_one[ 2 * i ] == 1 &&
        ptr_one[ 2 * i + 1 ] != 1 ) {
      rm = false;
      for( int j = 0; j < 6; ++ j ) {
        nidx = ptr_nidx[ 6 * i + j ];
        if( nidx != NA_INTEGER ) {
          if( ptr_one[ 2 * ( nidx - 1 ) ] == 0 ) {
            rm = true;
            break;
          }
        } else {
          rm = true;
          break;
        }
      }
      if( rm ) {
        outer_old.push_back( i + 1 );
        ptr_one[ 2 * i + 1 ] = 1;
      }
    }
  }
  while( outer_old.size() != 0 ) {
    // Peel outer part
    for( list<int>::const_iterator it = outer_old.begin();
         it != outer_old.end(); ++ it ) {
      idx = *it;
      ptr_one[ 2 * ( idx - 1 ) ] = 0;
      last = idx;
    }
    // Find new outer part
    for( list<int>::const_iterator it = outer_old.begin();
         it != outer_old.end(); ++ it ) {
      idx = *it;
      for( int i = 0; i < 6; ++ i ) {
        nidx = ptr_nidx[ 6 * ( idx - 1 ) + i ];
        if( nidx != NA_INTEGER ) {
          if( ptr_one[ 2 * ( nidx - 1 ) ] == 1 && 
              ptr_one[ 2 * nidx - 1 ] != 1 ) {
            outer_new.push_back( nidx );
            ptr_one[ 2 * nidx - 1 ] = 1; 
          }
        }
      }
    }
    outer_old = outer_new;
    clearVector( outer_new );
  }
  zeroVector( ptr_one, len );
  return last;
}