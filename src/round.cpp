#include "Rinternals.h"

#include "round.h"
#include "spread.h"
#include "radius.h"

double round( const int *ptr_seg1, const int &label1, 
              const int *ptr_seg2, const int &label2,
              const int &len, const int *ptr_nidx, 
              const int *ptr_aidx ) {
  vector<int> region;
  region.reserve( len );
  double roundness = 0;
  int idx = 0, nidx = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_seg2[ 2 * i ] == label2 ) {
      for( int j = 0; j < 6; ++ j ) {
        nidx = ptr_nidx[ 6 * i + j ];
        if( nidx != NA_INTEGER ) {
          if( ptr_seg2[ 2 * ( nidx - 1 ) ] != label2 &&
              ptr_seg1[ 2 * ( nidx - 1 ) ] == label1 ) {
            region.push_back( i + 1 );
            break;
          }
        }
      }
    }
  }
  double spread_idx = spread( region, ptr_aidx );
  double r2 = radius2D( region.size() );
  double r = radius( region.size() );
  
  roundness = spread_idx * r / r2;
  return roundness;
}