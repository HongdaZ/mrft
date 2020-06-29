#ifndef ADDREGION_H
#define ADDREGION_H

#include <list>
#include <vector>
#include <map>
#include <set>

using std::list;
using std::vector;
using std::map;
using std::set;

// add new tumor regions
void addRegion( const list<vector<double>> &region_parm, 
                const list<map<int,int>> &regions, 
                set<int> &tumor_labels,
                map<int, set<int>> &tumor_regions, 
                map<int, vector<double>> &tumor_parm );

#endif