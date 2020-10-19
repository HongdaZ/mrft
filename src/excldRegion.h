#ifndef EXCLDREGION_H
#define EXCLDREGION_H

#include <vector>

using std::vector;

// Remove the region from ptr_seg1 if there is no adjacent voxel with 
// ptr_seg2 == label
void excldRegion( const vector<int> &region, const int *ptr_nidx,
                  int *ptr_seg1,
                  const int *ptr_seg2, const int &label );

#endif