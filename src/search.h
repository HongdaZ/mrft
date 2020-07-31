#ifndef SEARCH_H
#define SEARCH_H

#include <list> 
#include <vector>

using std::list;
using std::vector;

void search( const int n_region, const bool init, int *label,
             const int *nidx, vector<int>& front, const int & l,
             vector<int>& region, list<int> &tumor_nbr, 
             bool &early_return );

#endif
