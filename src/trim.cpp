#include "Rinternals.h"

#include "trim.h"
#include "peel.h"
#include "pad2zero.h"
#include "cnctRegion.h"
#include "zeroVector.h"
#include "grow.h"

void trim( int *ptr_tumor, const int *ptr_nidx, 
           const int *ptr_aidx, vector<int> &region,
           const int &len, const double &s_trim ) {
  int *ptr_whole = new int[ 2 * len ]();
  int *ptr_res = new int[ 2 * len ]();
  int *ptr_one = new int[ 2 * len ]();
  int *ptr_keep = new int[ 2 * len ]();
  int *ptr_remain = new int[ 2 * len ]();
  int last = 0, n_tumor = 0, idx = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_tumor[ 2 * i ] == 1 ) {
      ptr_whole[ 2 * i ] = ptr_tumor[ 2 * i ];
      ++ n_tumor;
    }
  }
  // Let peel apples
  while( n_tumor != 0 ) {
    for( int i = 0; i < len; ++ i ) {
      if( cnctRegion( i + 1, ptr_nidx, ptr_whole, 
                      ptr_whole, 1, region ) ) {
        pad2zero( ptr_whole, region );
        zeroVector( ptr_one, len );
        for( vector<int>::const_iterator it = region.begin();
             it != region.end(); ++ it ) {
          idx = *it;
          ptr_one[ 2 * ( idx - 1 ) ] = 1;
        }
        last = peel( ptr_one, len, ptr_nidx );
        for( vector<int>::const_iterator it = region.begin();
             it != region.end(); ++ it ) {
          idx = *it;
          ptr_one[ 2 * ( idx - 1 ) ] = 1;
        }
        
        int n_ptr_one = 0;
        for( int j = 0; j < len; ++ j ) {
          if( ptr_one[ 2 * j ] == 1 ) {
            ++ n_ptr_one;
          }
        }
        grow( last, n_tumor, len, region, ptr_nidx, 
              ptr_aidx, ptr_whole, ptr_res, ptr_one,
              ptr_keep, ptr_remain, s_trim );
      }
    }
  }
  for( int i = 0; i < len; ++ i ) {
    if( ptr_res[ 2 * i ] == 0 ) {
      ptr_tumor[ 2 * i ] = 0;
    }
  }
  delete [] ptr_whole;
  delete [] ptr_res;
  delete [] ptr_one;
  delete [] ptr_keep;
  delete [] ptr_remain;
}