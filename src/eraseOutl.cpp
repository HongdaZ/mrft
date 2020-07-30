#include "eraseOutl.h"

// erase single outlier
void eraseOutl( const int out_label, 
                vector<int> &outl_labels,
                map<int, vector<double>> &outl_parm,
                int &n_outl ) {
  outl_labels[ out_label - 1 ] = 0;
  outl_parm.erase( out_label );
  -- n_outl;
  return;
}