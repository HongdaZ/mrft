#include "newTumorLabel.h"

// find a new tumor label for current voxel
int newTumorLabel( const int order, const set<int> &tumor_labels ) {
  int count = 0;
  int region_label = - 3;
  while( count < order ) {
    -- region_label;
    if( tumor_labels.find( region_label ) == tumor_labels.end() ) {
      ++ count;
    }
  }
  return region_label;
}