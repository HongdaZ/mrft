#ifndef SEARCH_H
#define SEARCH_H

#include <Rinternals.h>

#include <list> 
#include <queue>

using std::list;
using std::queue;

extern "C" {
  void search( const int n_region, const bool init, int *label,
               const int *nidx, queue<int>& front, const int & l,
               list<int>& region, list<int> &tumor_nbr, 
               bool &early_return );

} // extern "C"
#endif
