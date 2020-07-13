#include "addOutl.h"
// add new outlier
void addOutl( const int new_out_label,
              const vector<double> new_out_parm, 
              set<int> &outl_labels,
              map<int, vector<double>> &outl_parm ) {
  
  outl_labels.insert( new_out_label );
  outl_parm[ new_out_label ] = new_out_parm;
  
  return;
}