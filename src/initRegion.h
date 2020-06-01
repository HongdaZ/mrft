#ifndef INITREGION_H
#define INITREGION_H
#include <list>
#include <map>

using std::list;
using std::map;

map<int, list<int>> initRegion( int *ptr_seg, const int *ptr_nidx, int len, 
                                list<int> &labels );

#endif