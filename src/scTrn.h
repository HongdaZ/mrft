#ifndef SCTRN_H
#define SCTRN_H

#include <list>
#include <map>

using std::list;
using std::map;

int scTrn( list<int> &labels, list<list<int>> &regions, 
           const list<int> &tumor_labels,
           map<int, list<int>> &tumor_regions, 
           int *label, const int *nidx, int start );

#endif
