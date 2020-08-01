#ifndef GETREGION_H
#define GETREGION_H

#include <vector>

using std::vector;
// get the region with label = curr_label
void getRegion( vector<int> &region, const int &len_region, 
                const int &label, const int *ptr_seg, const int &len );

#endif