#ifndef SEARCH_H
#define SEARCH_H
#include <Rinternals.h>
#include <stack> 
std::stack<int> search( SEXP label, SEXP nidx, 
                       const int start, const double l );
#endif