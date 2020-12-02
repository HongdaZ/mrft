#include "Rinternals.h"

#include "trim.h"
#include "pad2zero.h"
#include "cnctRegion.h"
#include "excldRegion.h"
#include "zeroVector.h"
#include "spread.h"
#include "pTNbr.h"
#include "clearVector.h"
#include "radius.h"
#include "peel.h"

void trim( int *ptr_tumor, const int *ptr_nidx, vector<int> &region,
           const int &len ) {
  int *ptr_whole = new int[ 2 * len ]();
  int *ptr_remain = new int[ 2 * len ]();
  int last = 0;
  for( int i = 0; i < len; ++ i ) {
    ptr_whole[ 2 * i ] = ptr_tumor[ 2 * i ];
  }
  // last = index of last peeled voxel
  last = peel( ptr_whole, len, ptr_nidx );
  
  delete [] ptr_whole;
  delete [] ptr_remain;
}