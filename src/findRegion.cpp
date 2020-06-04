#include <R.h>
#include <Rinternals.h>
#include <list>
#include <queue>
#include "search.h"
#include "findRegion.h"

set<int> findRegion( int *ptr_label, const int *ptr_nidx , int start ) {

  queue<int> front;
  set<int> region;
  
  int l = ptr_label[ 2 * ( start - 1 ) ];
  front.push( start );
  region.insert( start );
  ptr_label[ 2 * start - 1 ] = 1;
  search( ptr_label, ptr_nidx, front, l, region );
  for( set<int>::iterator it = region.begin(); it != region.end(); ++ it ) {
    ptr_label[ 2 * *it - 1 ] = 0;
  }

  return region;
}