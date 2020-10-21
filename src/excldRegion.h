#ifndef EXCLDREGION_H
#define EXCLDREGION_H

#include <vector>

using std::vector;

// Remove the region from ptr_seg1 if there is no adjacent voxel with 
// ptr_seg2 == label
void excldRegion( const vector<int> &region, const int *ptr_nidx,
                  int *ptr_seg1,
                  const int *ptr_seg2, const int &label );
// Remove the region from ptr_seg1 if the number of voxels with 
// ptr_seg2 == label is less than size
void excldRegion( const vector<int> &region,
                  int *ptr_seg1,
                  const int *ptr_seg2, 
                  const int &label, const int &size );
// Remove tumor regions with size < size.
void excldRegion( const vector<int> &region,
                  int *ptr_seg, const int &size );

#endif