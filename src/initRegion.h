#ifndef INITREGION_H
#define INITREGION_H
#include <list>
#include <map>

using std::list;
using std::map;

void initRegion( int *ptr_seg, const int *ptr_nidx, int len, 
                 map<int, list<int>> &tumor, list<int> &labels );

#endif