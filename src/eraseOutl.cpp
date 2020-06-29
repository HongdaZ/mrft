#include "eraseOutl.h"

void eraseOutl( const int out_label, 
                set<int> &outl_labels,
                map<int, vector<double>> &outl_parm ) {
  outl_labels.erase( out_label );
  outl_parm.erase( out_label );
  return;
}