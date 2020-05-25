#ifndef SEARCH_H
#define SEARCH_H
#include <Rinternals.h>
#include <stack> 
#include <queue>
using std::stack;
using std::queue;
extern "C" {
  void search( double *label, int *nidx, 
               queue<int>& front, const int & l, stack<int>& region );
} // extern "C"
#endif