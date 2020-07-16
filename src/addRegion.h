#ifndef ADDREGION_H
#define ADDREGION_H

#include <list>
#include <vector>
#include <map>

using std::list;
using std::vector;
using std::map;

// add new tumor regions
void addRegion( int *ptr_seg, 
                const list<vector<double>> &region_parm, 
                const list<list<int>> &regions, 
                list<int> &tumor_labels,
                map<int, list<int>> &tumor_regions, 
                map<int, vector<double>> &tumor_parm );

#endif