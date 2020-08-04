#ifndef FINDREGION_H
#define FINDREGION_H

#include <vector>
#include <list>

using std::list;
using std::vector;

void findRegion( vector<int> &region, vector<int> &front,
                 const int &n_region, int *ptr_label, const int *ptr_nidx, 
                 const bool &init, list<int> &tumor_nbr, bool &early_return, 
                 const int &start );

#endif