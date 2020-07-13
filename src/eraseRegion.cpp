#include "eraseRegion.h"

// erase single tumor region
void eraseRegion( const int tumor_label, 
                  set<int> &tumor_labels,
                  map<int, set<int>> &tumor_regions, 
                  map<int, vector<double>> &tumor_parm ) {
  tumor_labels.erase( tumor_label );
  tumor_regions.erase( tumor_label );
  tumor_parm.erase( tumor_label );
  return;
}