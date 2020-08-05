#ifndef GETREGION_H
#define GETREGION_H

#include <vector>
#include <list>

using std::vector;
using std::list;

// get the region with label = label
void getRegion( vector<int> &region, const int &label, const int *ptr_seg,
                const int &len );
void getRegion( vector<int> &region, const list<int> &t_region );

#endif