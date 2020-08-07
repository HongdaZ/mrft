#ifndef GETREGION_H
#define GETREGION_H

#include <Rinternals.h>

#include <vector>

using std::vector;

// get region and region_label
void getRegion( int &region_label, vector<int> &region, 
                const vector<int> &regions_whole, 
                const vector<int> &regions_sub, 
                const int &order );
#endif