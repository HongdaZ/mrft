#include "grow.h"
#include "pad2zero.h"
#include "cnctRegion.h"
#include "zeroVector.h"
#include "spread.h"
#include "pTNbr.h"
#include "clearVector.h"
#include "radius.h"
#include "peel.h"
#include "nTumor.h"
#include "nSurface.h"

void grow( const int &last, int &n_tumor, const int &len, 
           vector<int> &region,
           const int *ptr_nidx, int *ptr_whole, int *ptr_res,
           int *ptr_one, int *ptr_keep, int *ptr_remain ) {
  list<int> outer;
  zeroVector( ptr_keep, len );
  zeroVector( ptr_remain, len );
  for( int i = 0; i < len; ++ i ) {
    ptr_remain[ 2 * i ] = ptr_one[ 2 * i ];
  }
  
}