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

void initRegion( vector<int> &region, vector<int> &front, int *ptr_seg, 
                 const int *ptr_nidx, int len, vector<int> &n_voxel, 
                 vector<int> &labels, int &n_tumor );

#endif