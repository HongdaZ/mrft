#ifndef ERASEOUTL_H
#define ERASEOUTL_H

#include <list>
#include <map>
#include <vector>

using std::list;
using std::map;
using std::vector;

// erase single outlier
void eraseOutl( const int out_label, 
                list<int> &outl_labels,
                map<int, vector<double>> &outl_parm );

#endif