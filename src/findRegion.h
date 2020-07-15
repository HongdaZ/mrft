#ifndef FINDREGION_H
#define FINDREGION_H
#include <set>
using std::set;

set<int> findRegion( int *ptr_label, const int *ptr_nidx, const bool init,
                     set<int> &tumor_nbr, bool &early_return,int start );

#endif