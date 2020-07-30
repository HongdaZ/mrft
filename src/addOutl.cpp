#include "addOutl.h"
// add new outlier
void addOutl( const int new_out_label,
              const vector<double> &new_out_parm, 
              vector<int> &outl_labels,
              map<int, vector<double>> &outl_parm,
              int &n_outl ) {
  
  outl_labels[ new_out_label - 1 ] = 1;
  outl_parm[ new_out_label ] = new_out_parm;
  ++ n_outl;
  return;
}