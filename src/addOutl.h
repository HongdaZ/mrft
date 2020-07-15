#ifndef ADDOUTL_H
#define ADDOUTL_H

#include <vector>
#include <map>
#include <set>

using std::vector;
using std::map;
using std::set;

// add new outlier
void addOutl( const int new_out_label,
              const vector<double> &new_out_parm, 
              set<int> &outl_labels,
              map<int, vector<double>> &outl_parm );

#endif