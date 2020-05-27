#ifndef SEARCH_H
#define SEARCH_H
#include <stack> 
#include <queue>
#include <Rinternals.h>
using std::stack;
using std::queue;

extern "C" {
  void search( int *label, int *nidx, 
               queue<int>& front, const int & l, stack<int>& region );

} // extern "C"
#endif
