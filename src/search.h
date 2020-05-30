#ifndef SEARCH_H
#define SEARCH_H

#include <Rinternals.h>

#include <list> 
#include <queue>

using std::list;
using std::queue;

extern "C" {
  void search( int *label, int *nidx, 
               queue<int>& front, const int & l, list<int>& region );

} // extern "C"
#endif
