#include "newTumorLabel.h"
#include "algorithm"

// find a new tumor label for current voxel
int newTumorLabel( const int order, const vector<int> &tumor_labels ) {
  int count = 0;
  int region_label;
  int i = 0;
  for( ; i < tumor_labels.size(); ++ i ) {
    if( tumor_labels[ i ] == 0 ) {
      ++ count;
      if( count == order ) {
        break;
      }
    }
  }
  region_label = - i - 4;
  return region_label;
}