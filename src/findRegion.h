#ifndef FINDREGION_H
#define FINDREGION_H
#include <list>
using std::list;

list<int> findRegion( const int n_region, int *ptr_label, 
                       const int *ptr_nidx, const bool init,
                       list<int> &tumor_nbr, bool &early_return,int start );

#endif