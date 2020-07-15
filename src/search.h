#ifndef SEARCH_H
#define SEARCH_H

#include <Rinternals.h>

#include <set> 
#include <queue>

using std::set;
using std::queue;

extern "C" {
  void search( const bool init, int *label, const int *nidx, 
               queue<int>& front, const int & l, set<int>& region,
               set<int> &tumor_nbr, bool &early_return );

} // extern "C"
#endif
