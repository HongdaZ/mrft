#include <R.h>
#include <Rinternals.h>
#include <list>
#include <queue>
#include "search.h"
#include "findRegion.h"

set<int> findRegion( const int n_region, int *ptr_label, const int *ptr_nidx, const bool init, 
                     set<int> &tumor_nbr, bool &early_return,
                     int start ) {

  queue<int> front;
  set<int> region;
  
  int l = ptr_label[ 2 * ( start - 1 ) ];
  front.push( start );
  region.insert( start );
  ptr_label[ 2 * start - 1 ] = 1;
  search( n_region, init, ptr_label, ptr_nidx, front, l, region, tumor_nbr,
          early_return );
  for( set<int>::iterator it = region.begin(); it != region.end(); ++ it ) {
    ptr_label[ 2 * *it - 1 ] = 0;
  }

  return region;
}