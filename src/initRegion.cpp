#include <R.h>
#include <Rinternals.h>

#include <map>
#include <list>

#include "findRegion.h"
#include "initRegion.h"

using::std::map;
using::std::list;

map<int, list<int>> initRegion( int *ptr_seg, const int *ptr_nidx, int len,
                                list<int> &labels ) {
  
  map<int, list<int>> tumor;
  list<int> region;
  int tumor_label = - 4;
  
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg[ 2 * i ] > 0 ) {
      region = findRegion( ptr_seg, ptr_nidx, i + 1 );
      tumor[ tumor_label ] = region;
      for( list<int>::iterator it = region.begin(); it != region.end(); 
      ++ it ) {
        ptr_seg[ 2 * ( *it - 1 ) ] = tumor_label;
      }
      labels.push_back( tumor_label );
      -- tumor_label;
    }
  }
  return tumor;
}
