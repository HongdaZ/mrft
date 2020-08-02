#ifndef GETREGION_H
#define GETREGION_H

#include <vector>

using std::vector;
// get the region with label = label
void getRegion( vector<int> &region, const int &label, const int *ptr_seg,
                const int &len );
void getRegion( vector<int> &region, const int &label, 
                const vector<int> &regions,
                const int &len, const int &row );

#endif