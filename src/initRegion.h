#ifndef INITREGION_H
#define INITREGION_H
#include <set>
#include <map>
#include <list>
#include <vector>

using std::set;
using std::map;
using std::list;
using std::vector;

void initRegion( int *ptr_seg, const int *ptr_nidx, int len, 
                 map<int, list<int>> &tumor, vector<int> &labels );

#endif