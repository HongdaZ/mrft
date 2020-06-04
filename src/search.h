#ifndef SEARCH_H
#define SEARCH_H

#include <Rinternals.h>

#include <set> 
#include <queue>

using std::set;
using std::queue;

extern "C" {
  void search( int *label, const int *nidx, 
               queue<int>& front, const int & l, set<int>& region );

} // extern "C"
#endif
