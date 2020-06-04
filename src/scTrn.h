#ifndef SCTRN_H
#define SCTRN_H

#include <list>
#include <map>

using std::list;
using std::map;

int scTrn( list< map<int, int>> &regions, 
           int *label, const int *nidx, int start );

#endif
