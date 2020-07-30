#ifndef ADDOUTL_H
#define ADDOUTL_H

#include <vector>
#include <map>

using std::vector;
using std::map;

// add new outlier
void addOutl( const int new_out_label,
              const vector<double> &new_out_parm, 
              vector<int> &outl_labels,
              map<int, vector<double>> &outl_parm,
              int &n_outl );

#endif