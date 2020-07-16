#ifndef ERASEREGION_H
#define ERASEREGION_H

#include <list>
#include <map>
#include <vector>

using std::list;
using std::map;
using std::vector;

// erase single tumor region
void eraseRegion( const int tumor_label, 
                 list<int> &tumor_labels,
                 map<int, list<int>> &tumor_regions, 
                 map<int, vector<double>> &tumor_parm );

#endif