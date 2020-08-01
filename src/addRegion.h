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
                const vector<double> &region_parm,
                const vector<double> &regions, const int &region_row,
                vector<int> &tumor_labels,
                vector<double> &tumor_parm,
                int &n_tumor );

#endif