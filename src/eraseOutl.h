#ifndef ERASEOUTL_H
#define ERASEOUTL_H

#include <map>
#include <vector>

using std::map;
using std::vector;

// erase single outlier
void eraseOutl( const int out_label, 
                vector<int> &outl_labels,
                map<int, vector<double>> &outl_parm,
                int &n_outl );

#endif