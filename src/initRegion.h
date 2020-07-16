#ifndef INITREGION_H
#define INITREGION_H
#include <set>
#include <map>
#include <list>

using std::set;
using std::map;
using std::list;

void initRegion( int *ptr_seg, const int *ptr_nidx, int len, 
                 map<int, list<int>> &tumor, list<int> &labels );

#endif