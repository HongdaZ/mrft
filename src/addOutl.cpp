#include "addOutl.h"
// add new outlier
void addOutl( const int new_out_label,
              const vector<double> &new_out_parm, 
              vector<int> &outl_labels,
              vector<double> &outl_parm,
              int &n_outl ) {
  int col = new_out_label - 1;
  outl_labels[ col ] = 1;
  for( int i = 0; i < 2; ++ i ) {
    outl_parm[ 2 * col + i ] = new_out_parm[ i ];
  }
  ++ n_outl;
  return;
}