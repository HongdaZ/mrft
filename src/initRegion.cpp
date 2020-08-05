#include <R.h>
#include <Rinternals.h>

#include <list>

#include "findRegion.h"
#include "initRegion.h"
#include "label2col.h"

void initRegion( vector<int> &region, vector<int> &front, 
                 list<list<int>> &tumor_regions, 
                 list<int> &tumor_regions_label, int *ptr_seg, 
                 const int *ptr_nidx, const int &len,
                 vector<int> &labels, int &n_tumor ) {
  int tumor_label = - 4;
  int n_region = 1;
  int idxcol;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg[ 2 * i ] > 0 ) {
      list<int> tumor_nbr;
      // a new tumor region
      list<int> new_region;
      tumor_regions_label.push_back( tumor_label );
      bool early_return;
      findRegion( region, front, n_region, ptr_seg, ptr_nidx, true, 
                  tumor_nbr, early_return, i + 1 );
      for( vector<int>::iterator it = region.begin(); it != region.end(); 
      ++ it ) {
        ptr_seg[ 2 * ( *it - 1 ) ] = tumor_label;
        new_region.push_back( *it );
      }
      tumor_regions.push_back( new_region );
      idxcol = label2col( tumor_label );
      labels[ idxcol ] = 1;
      -- tumor_label;
      ++ n_tumor;
    }
  }
  return;
}
