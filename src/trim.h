#ifndef TRIM_H
#define TRIM_H

#include <vector>

using std::vector;

// Remove narrow path
void trim( int *ptr_tumor, const int *ptr_nidx, 
           const int *ptr_aidx, vector<int> &region,
           const int &len, const double &s_trim );
#endif