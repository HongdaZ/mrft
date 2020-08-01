#include "eraseRegion.h"

// erase single tumor region
void eraseRegion( const int tumor_label, 
                  vector<int> &tumor_labels,
                  vector<double> &tumor_parm,
                  int &n_tumor ) {
  int idxcol = - tumor_label - 4;
  tumor_labels[ idxcol ] = 0;
  for( int i = 0; i < 8; ++ i ) {
    tumor_parm[ 8 * idxcol + i ] = 0;
  }
  -- n_tumor;
  return;
}