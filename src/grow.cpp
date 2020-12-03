#include "Rinternals.h"

#include "grow.h"
#include "pad2zero.h"
#include "cnctRegion.h"
#include "zeroVector.h"
#include "spread.h"
#include "pTNbr.h"
#include "clearVector.h"
#include "radius.h"
#include "nTumor.h"
#include "nSurface.h"

void grow( const int &last, int &n_tumor, const int &len, 
           vector<int> &region,
           const int *ptr_nidx, const int *ptr_aidx,
           int *ptr_whole, int *ptr_res,
           int *ptr_one, int *ptr_keep, int *ptr_remain ) {
  list<int> outer_old, outer_new;
  int idx = 0, nidx = 0, spread_keep = 0, spread_remain = 0;
  int n_surface = 0, n_tnbr = 0;
  double p_tnbr = 0, r = 0;
  bool keep;
  zeroVector( ptr_keep, len );
  zeroVector( ptr_remain, len );
  for( int i = 0; i < len; ++ i ) {
    ptr_remain[ 2 * i ] = ptr_one[ 2 * i ];
  }
  outer_old.push_back( last );
  keep = true;
  while( outer_old.size() > 0 ) {
    // Add outer part to keep
    for( list<int>::const_iterator it = outer_old.begin();
         it != outer_old.end(); ++ it ) {
      idx = *it;
      ptr_keep[ 2 * ( idx - 1 ) ] = 1;
      ptr_remain[ 2 * ( idx - 1 ) ] = 0;
      -- n_tumor;
    }
    cnctRegion( idx, ptr_nidx, ptr_keep, ptr_keep, 1, region );
    pad2zero( ptr_keep, region );
    spread_keep = spread( region, ptr_aidx );
    if( spread_keep > 4 ) {
      keep = false;
      break;
    }
    n_surface = nSurface( ptr_keep, len, ptr_nidx );
    for( int i = 0; i < len; ++ i ) {
      if( cnctRegion( i + 1, ptr_nidx, ptr_remain, ptr_remain, 1,
                      region ) ) {
        p_tnbr = pTNbr( region, ptr_keep, 1, ptr_nidx );
        n_tnbr = p_tnbr * region.size();
        spread_remain = spread( region, ptr_aidx );
        r = radius( region.size() );
        if( n_tnbr < ( double )n_surface / 54 &&
            1 / p_tnbr > 4 * r &&
            spread_remain > 4 ) {
          // Remove from ptr_one and ptr_remain
          for( vector<int>::const_iterator it = region.begin();
               it != region.end(); ++ it ) {
            idx = *it;
            if( idx != 0 ) {
              ptr_one[ 2 * ( idx - 1 ) ] = 0;
              ptr_remain[ 2 * ( idx - 1 ) ] = 0;
            }
          }
        }
      }
    }
    pad2zero( ptr_remain, len );
    // Find new outer part
    for( list<int>::const_iterator it = outer_old.begin();
         it != outer_old.end(); ++ it ) {
      idx = *it;
      for( int i = 0; i < 6; ++ i ) {
        nidx = ptr_nidx[ 6 * ( idx - 1 ) + i ];
        if( nidx != NA_INTEGER ) {
          if( ptr_remain[ 2 * ( nidx - 1 ) ] == 1 ) {
            outer_new.push_back( nidx );
          }
        }
      }
    }
    outer_old = outer_new;
    clearVector( outer_new );
    
  }
  // Keep the tumor region
  if( keep ) {
    for( int i = 0; i < len; ++ i ) {
      if( ptr_keep[ 2 * i ] == 1 ) {
        ptr_res[ 2 * i ] = 1;
      }
    }
  }
  for( int i = 0; i < len; ++ i ) {
    if( ptr_keep[ 2 * i ] == 1 ) {
      ptr_whole[ 2 * i ] = 0;
    }
  }
}