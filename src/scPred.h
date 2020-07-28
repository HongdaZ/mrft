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
            const list<int> &tumor_labels,
            map<int, list<int>> &tumor_regions,
            int *ptr_label, const int *ptr_nidx, int start );

#endif