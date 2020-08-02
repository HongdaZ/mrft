#include "addRegion.h"
// add new tumor regions
void addRegion( int *ptr_seg, 
                const vector<double> &region_parm,
                const vector<int> &regions, const int &region_row,
                vector<int> &tumor_labels,
                vector<double> &tumor_parm,
                int &n_tumor ) {
  int nrow = 1 + 2 + 6;
  int ncol = region_parm.size() / nrow;
  int new_label;
  int idxcol;
  for( int i = 0; i < ncol; ++ i ) {
    new_label = region_parm[ nrow * i ];
    idxcol = - new_label - 4;
    tumor_labels[ idxcol ] = 1;
    for( int j = 0; j < 8; ++ j ) {
      tumor_parm[ 8 * idxcol + j ] = region_parm[ nrow * i + j + 1 ];
    }
    ++ n_tumor;
  }
  for( int i = 1; i <= regions.size() / 2; ++ i ) {
    ptr_seg[ 2 * ( i - 1 ) ] = regions[ 2 * ( i - 1 ) + region_row ];
  }
  return;
}