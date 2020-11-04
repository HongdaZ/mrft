#ifndef SPREAD_H
#define SPREAD_H

#include <vector>

using std::vector;

// Measure the spread of the voxels in region
double spread( const vector<int> &region, const int *ptr_aidx );

#endif