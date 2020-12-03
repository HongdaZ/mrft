#include "Rinternals.h"

#include "trim.h"
#include "peel.h"
#include "pad2zero.h"

void trim( int *ptr_tumor, const int *ptr_nidx, vector<int> &region,
           const int &len ) {
  int *ptr_whole = new int[ 2 * len ]();
  int *ptr_res = new int[ 2 * len ]();
  int *ptr_one = new int[ 2 * len ]();
  int *ptr_keep = new int[ 2 * len ]();
  int *ptr_remain = new int[ 2 * len ]();
  int last = 0;
  for( int i = 0; i < len; ++ i ) {
    ptr_whole[ 2 * i ] = ptr_tumor[ 2 * i ];
  }
  // last = index of last peeled voxel
  
  delete [] ptr_whole;
  delete [] ptr_res;
  delete [] ptr_one;
  delete [] ptr_keep;
  delete [] ptr_remain;
}