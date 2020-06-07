#include <R.h>
#include <Rinternals.h>

#include <list>

#include "findRegion.h"
#include "initRegion.h"

using std::map;
using std::set;
using std::list;

void initRegion( int *ptr_seg, const int *ptr_nidx, int len,
                                map<int, set<int>> &tumor,
                                set<int> &labels ) {
  
  set<int> region;
  int tumor_label = - 4;
  
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg[ 2 * i ] > 0 ) {
      region = findRegion( ptr_seg, ptr_nidx, i + 1 );
      tumor[ tumor_label ] = region;
      // // Debug
      // if( tumor_label < -20 ) {
      //   for( set<int>::iterator it = region.begin(); it != region.end(); ++ it ) {
      //     Rprintf("region index: %d\n", *it );
      //   }
      //   for( set<int>::iterator it = tumor[ tumor_label ].begin(); it != tumor[ tumor_label ].end(); ++ it ) {
      //     Rprintf("tumor regions index: %d\n", *it );
      //   }
      // }
      // //////////////////////////////////////////////////////////
      for( set<int>::iterator it = region.begin(); it != region.end(); 
      ++ it ) {
        ptr_seg[ 2 * ( *it - 1 ) ] = tumor_label;
      }
      labels.insert( tumor_label );
      -- tumor_label;
    }
  }
  return;
}
