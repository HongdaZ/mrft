#ifndef SCTRN_H
#define SCTRN_H

#include <list>
#include <map>
#include <vector>

using std::list;
using std::map;
using std::vector;

int scTrn( list<int> &labels, list<list<int>> &regions, 
           const vector<int> &tumor_labels,
           map<int, list<int>> &tumor_regions, 
           int *label, const int *nidx, int start );

#endif
