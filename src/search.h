#ifndef SEARCH_H
#define SEARCH_H
#include <Rinternals.h>
#include <stack> 
#include <queue>
using std::stack;
using std::queue;

void search( int *label, int *nidx, 
             queue<int>& front, const int & l, stack<int>& region );
#endif