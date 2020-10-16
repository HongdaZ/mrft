#ifndef REGIONS_H
#define REGIONS_H

#include <vector>
#include <list>

using std::vector;
using std::list;

list<vector<int>> regions( int *ptr_seg, const int &label,
                           const int *ptr_nidx,
                           const int *ptr_aidx );

#endif