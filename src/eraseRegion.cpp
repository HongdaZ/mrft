#include "eraseRegion.h"

// erase single tumor region
void eraseRegion( const int tumor_label, 
                  vector<int> &tumor_labels,
                  map<int, list<int>> &tumor_regions, 
                  map<int, vector<double>> &tumor_parm,
                  int &n_tumor ) {
  tumor_labels[ - tumor_label - 4 ] = 0;
  tumor_regions.erase( tumor_label );
  tumor_parm.erase( tumor_label );
  -- n_tumor;
  return;
}