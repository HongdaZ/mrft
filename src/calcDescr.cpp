#include "calcDescr.h"
#include "cHull.h"
#include "polygonArea.h"
#include "spread.h"
#include "roundness.h"

void calcDescr( const int &size_region,
                int *ptr_seg_copy, const int &len, 
                vector<int> &region,
                const vector<list<vector<int>>> &slices, 
                vector<double> &solidity, 
                vector<double> &avg_spread, 
                vector<double> &avg_round,
                const int *ptr_aidx, const int *ptr_nidx ) {
  // Axial, sagittal and coronal
  vector<int> row{ 2, 3, 1, 3, 1, 2 };
  double hull_area;
  int len_slice;
  for( int i = 0; i < 3; ++ i ) {
    solidity[ i ] = 0;
    avg_spread[ i ] = 0;
    avg_round[ i ] = 0;
    const list<vector<int>> &p_slice = slices[ i ];
    for( list<vector<int>>::const_iterator it_s = p_slice.begin();
         it_s != p_slice.end(); ++ it_s ) {
      
      len_slice = it_s->size() / 4;
      vector<int> chull_input( 2 * len_slice );
      for( int j = 0; j < len_slice; ++ j ) {
        chull_input[ 2 * j ] = ( *it_s )[ 4 * j + row[ 2 * i ] ];
        chull_input[ 2 * j + 1 ] = 
          ( *it_s )[ 4 * j + row[ 2 * i + 1 ] ];
      }
      vector<int> hull = cHull( chull_input );
      hull_area = polygonArea( hull );
      solidity[ i ] += hull_area;
      avg_spread[ i ] += spread( *it_s, ptr_seg_copy,
                            len, ptr_nidx, i );
      avg_round[ i ] += roundness( *it_s, i, ptr_aidx );

    }
    solidity[ i ] /= size_region;
    avg_spread[ i ] /= p_slice.size();
    avg_round[ i ] /= p_slice.size();
  }
}