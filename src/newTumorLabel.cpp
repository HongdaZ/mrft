#include "newTumorLabel.h"
#include "algorithm"

using std::find;

// find a new tumor label for current voxel
int newTumorLabel( const int order, const list<int> &tumor_labels ) {
  int count = 0;
  int region_label = - 3;
  while( count < order ) {
    -- region_label;
    list<int>::const_iterator it = find( tumor_labels.begin(),
                                         tumor_labels.end(), region_label );
    if( it == tumor_labels.end() ) {
      ++ count;
    }
  }
  return region_label;
}