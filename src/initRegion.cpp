#include <R.h>
#include <Rinternals.h>

#include <list>

#include "findRegion.h"
#include "initRegion.h"

void initRegion( int *ptr_seg, const int *ptr_nidx, int len,
                 map<int, list<int>> &tumor, vector<int> &labels ) {
  
  list<int> region;
  int tumor_label = - 4;
  int n_region = 1;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg[ 2 * i ] > 0 ) {
      list<int> tumor_nbr;
      bool early_return;
      region = findRegion( n_region, ptr_seg, ptr_nidx, true, 
                           tumor_nbr, early_return, i + 1 );
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
      for( list<int>::iterator it = region.begin(); it != region.end(); 
      ++ it ) {
        ptr_seg[ 2 * ( *it - 1 ) ] = tumor_label;
      }
      labels[ - tumor_label - 4 ] = 1;
      -- tumor_label;
    }
  }
  return;
}
