#ifndef NEWTUMORLABEL_H
#define NEWTUMORLABEL_H

#include <set>

using std::set;

// find an outlier label for current voxel
int newTumorLabel( const int order, const set<int> &tumor_labels );

#endif