#include <R.h>
#include <Rinternals.h>

#include <vector>

#include "search.h"
#include "findRegion.h"

#include "clearVector.h"

using std::vector;

void findRegion( vector<int> &region, vector<int> &front,
                 const int n_region, int *ptr_label, const int *ptr_nidx, 
                 const bool init, list<int> &tumor_nbr, bool &early_return, 
                 int start ) {
  clearVector( front );
  clearVector( region );
  
  int l = ptr_label[ 2 * ( start - 1 ) ];
  front.push_back( start );
  region.push_back( start );
  ptr_label[ 2 * start - 1 ] = 1;
  search( n_region, init, ptr_label, ptr_nidx, front, l, region, tumor_nbr,
          early_return );
  for( vector<int>::iterator it = region.begin(); it != region.end(); ++ it ) {
    ptr_label[ 2 * *it - 1 ] = 0;
  }

  return;
}