#include "eraseOutl.h"

// erase single outlier
void eraseOutl( const int out_label, 
                list<int> &outl_labels,
                map<int, vector<double>> &outl_parm ) {
  outl_labels.remove( out_label );
  outl_parm.erase( out_label );
  return;
}