#ifndef NEWTUMORLABEL_H
#define NEWTUMORLABEL_H

#include <list>

using std::list;

// find an outlier label for current voxel
int newTumorLabel( const int order, const list<int> &tumor_labels );

#endif