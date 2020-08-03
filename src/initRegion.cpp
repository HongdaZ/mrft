#include <R.h>
#include <Rinternals.h>

#include <list>

#include "findRegion.h"
#include "initRegion.h"

void initRegion( vector<int> &region, vector<int> &front, int *ptr_seg, 
                 const int *ptr_nidx, int len, vector<int> &n_voxel, 
                 vector<int> &labels, int &n_tumor ) {
  int tumor_label = - 4;
  int n_region = 1;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg[ 2 * i ] > 0 ) {
      list<int> tumor_nbr;
      bool early_return;
      findRegion( region, front, n_region, ptr_seg, ptr_nidx, true, 
                  tumor_nbr, early_return, i + 1 );
      n_voxel[ - tumor_label - 4 ] = region.size();
      for( vector<int>::iterator it = region.begin(); it != region.end(); 
      ++ it ) {
        ptr_seg[ 2 * ( *it - 1 ) ] = tumor_label;
      }
      labels[ - tumor_label - 4 ] = 1;
      -- tumor_label;
      ++ n_tumor;
    }
  }
  return;
}
