#ifndef INITREGION_H
#define INITREGION_H
#include <set>
#include <map>

using std::set;
using std::map;

void initRegion( int *ptr_seg, const int *ptr_nidx, int len, 
                 map<int, set<int>> &tumor, set<int> &labels );

#endif