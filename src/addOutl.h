#ifndef ADDOUTL_H
#define ADDOUTL_H

#include <vector>
#include <map>
#include <list>

using std::vector;
using std::map;
using std::list;

// add new outlier
void addOutl( const int new_out_label,
              const vector<double> &new_out_parm, 
              list<int> &outl_labels,
              map<int, vector<double>> &outl_parm );

#endif