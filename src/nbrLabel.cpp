#include <R.h>
#include <Rinternals.h>

#include <vector>

#include "nbrLabel.h"

using std::vector;
// find the labels of the neighbors
vector<int> nbrLabel( const int curr_idx, 
                 const int *ptr_seg, const int *ptr_nidx ) {
  
  vector<int> nbr_label( 6 );
  int nbr_idx = 0;

  for( int i = 0; i < 6; ++ i ) {
    nbr_idx = ptr_nidx[ 6 * ( curr_idx - 1 ) + i ];
    if( nbr_idx != NA_INTEGER ) {
      nbr_label[ i ] = ptr_seg[ 2 * ( nbr_idx - 1 ) ];
    } else {
      nbr_label[ i ] = NA_INTEGER;
    }
  }
  return nbr_label;
}
