#ifndef FINDOUTLABEL_H
#define FINDOUTLABEL_H

#include <vector>

using std::vector;

// find an outlier label for current voxel
int findOutLabel( const int idx, const int *ptr_seg, 
                  vector<int> &outl_labels );

#endif