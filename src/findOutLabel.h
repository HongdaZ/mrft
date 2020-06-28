#ifndef FINDOUTLABEL_H
#define FINDOUTLABEL_H

#include <set>

using std::set;

// find an outlier label for current voxel
int findOutLabel( const int idx, const int *ptr_seg, 
                 const set<int> &outl_labels );

#endif