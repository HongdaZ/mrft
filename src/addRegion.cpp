#include "addRegion.h"
#include "label2col.h"

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
  int cidx;
  for( int i = 0; i < ncol; ++ i ) {
    new_label = region_parm[ nrow * i ];
    cidx = label2col( new_label );
    tumor_labels[ cidx ] = 1;
    for( int j = 0; j < 8; ++ j ) {
      tumor_parm[ 8 * cidx + j ] = region_parm[ nrow * i + j + 1 ];
    }
    ++ n_tumor;
  }
  for( int i = 1; i <= regions.size() / 2; ++ i ) {
    new_label = regions[ 2 * ( i - 1 ) + region_row ];
    if( new_label != 0 ) {
      ptr_seg[ 2 * ( i - 1 ) ] = new_label;
    }
    
  }
  return;
}
// add a new single voxel tumor region
void addRegion( int *ptr_seg, 
                const vector<double> &region_parm,
                const vector<int> region,
                vector<int> &tumor_labels,
                vector<double> &tumor_parm,
                int &n_tumor ) {
  int nrow = 1 + 2 + 6;
  int ncol = region_parm.size() / nrow;
  int new_label;
  int cidx;
  for( int i = 0; i < ncol; ++ i ) {
    new_label = region_parm[ nrow * i ];
    cidx = label2col( new_label );
    tumor_labels[ cidx ] = 1;
    for( int j = 0; j < 8; ++ j ) {
      tumor_parm[ 8 * cidx + j ] = region_parm[ nrow * i + j + 1 ];
    }
    ++ n_tumor;
  }
  ptr_seg[ 2 * ( region[ 0 ] - 1 ) ] = new_label;
  return;
}