#ifndef FINDOUTLABEL_H
#define FINDOUTLABEL_H

#include <list>

using std::list;

// find an outlier label for current voxel
int findOutLabel( const int idx, const int *ptr_seg, 
                 const list<int> &outl_labels );

#endif