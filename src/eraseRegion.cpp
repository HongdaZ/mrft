#include "eraseRegion.h"

// erase single tumor region
void eraseRegion( const int tumor_label, 
                  list<int> &tumor_labels,
                  map<int, list<int>> &tumor_regions, 
                  map<int, vector<double>> &tumor_parm ) {
  tumor_labels.remove( tumor_label );
  tumor_regions.erase( tumor_label );
  tumor_parm.erase( tumor_label );
  return;
}