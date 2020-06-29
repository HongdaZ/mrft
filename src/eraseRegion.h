#ifndef ERASEREGION_H
#define ERASEREGION_H

#include <set>
#include <map>
#include <vector>

using std::set;
using std::map;
using std::vector;

// find an outlier label for current voxel
void eraseRegion( const int tumor_label, 
                 set<int> &tumor_labels,
                 map<int, set<int>> &tumor_regions, 
                 map<int, vector<double>> &tumor_parm );

#endif