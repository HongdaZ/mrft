#ifndef SCTRN_H
#define SCTRN_H

#include <list>
#include <map>
#include <set>

using std::list;
using std::map;
using std::set;

int scTrn( list< map<int, int>> &regions, 
           const set<int> &tumor_labels,
           map<int, set<int>> &tumor_regions, 
           int *label, const int *nidx, int start );

#endif
