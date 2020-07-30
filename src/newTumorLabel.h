#ifndef NEWTUMORLABEL_H
#define NEWTUMORLABEL_H

#include <vector>

using std::vector;

// find an outlier label for current voxel
int newTumorLabel( const int order, const vector<int> &tumor_labels );

#endif