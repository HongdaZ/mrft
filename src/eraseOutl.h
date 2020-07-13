#ifndef ERASEOUTL_H
#define ERASEOUTL_H

#include <set>
#include <map>
#include <vector>

using std::set;
using std::map;
using std::vector;

// erase single outlier
void eraseOutl( const int out_label, 
                 set<int> &outl_labels,
                 map<int, vector<double>> &outl_parm );

#endif