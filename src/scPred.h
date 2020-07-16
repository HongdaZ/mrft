#ifndef SCPRED_H
#define SCPRED_H

#include <list>
#include <map>
#include <set>

using std::list;
using std::map;
using std::set;

// split or combine for prediction
int scPred( list<int> &labels, list<list<int>> &regions, 
           const set<int> &tumor_labels,
           map<int, set<int>> &tumor_regions, 
           int *label, const int *nidx, int start );

#endif