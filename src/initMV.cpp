#include <R.h>
#include <Rinternals.h>

#include "initMV.h"

// initialize yln_, yln_i and yl_
void initMV( const list<int> &region, double *yln_, double *yln_i,
             double *yl_, double &sum_y, const double *ptr_intst, 
             const int *ptr_nidx,
             const double *ptr_nintst, const int *ptr_seg, 
             const int curr_label ) {
  int nrow = region.size();
  int ncol = 6;
  list<int>::const_iterator it = region.begin();
  int idx;
  for( int j = 0; j < nrow; ++ it, ++ j ) {
    idx = *it;
    yl_[ j ] = ptr_intst[ idx - 1 ];
    sum_y += yl_[ j ];
  }
  for( int i = 0; i < 6; ++ i ) {
    it = region.begin();
    for( int j = 0; j < nrow; ++ it, ++ j ) {
      idx = *it;
      int nidx =  ptr_nidx[ 6 * ( idx - 1 ) + i ];
      if( nidx != NA_INTEGER ) {
        // healthy
        if( curr_label < 0 && curr_label > - 4 ) {
          int nlabel = ptr_seg[ 2 * ( nidx - 1 ) ];
          if( nlabel == curr_label ) {
            yln_[ i * nrow + j ] = ptr_nintst[ 6 * ( idx - 1 ) + i ];
            yln_i[ i * nrow + j ] = 1;
          } else {
            yln_[ i * nrow + j ] = 0;
            yln_i[ i * nrow + j ] = 0;
          }
          // tumor
        } else {
          int nlabel = ptr_seg[ 2 * nidx - 1 ];
          if( nlabel == 3 ) {
            yln_[ i * nrow + j ] = ptr_nintst[ 6 * ( idx - 1 ) + i ];
            yln_i[ i * nrow + j ] = 1;
          } else {
            yln_[ i * nrow + j ] = 0;
            yln_i[ i * nrow + j ] = 0;
          }
        } 
      } else {
        yln_[ i * nrow + j ] = 0;
        yln_i[ i * nrow + j ] = 0;
      }
    }
  }
  return;
}
// initialize yl and yln
void initMV( double *yl, const double *yl_, double *yln, const double *yln_,
             const double *yln_i, const int nrow, const double mu ) {
  int ncol = 6;
  for( int j = 0; j < nrow; ++ j ) {
    yl[ j ] = yl_[ j ] - mu;
  }
  for( int i = 0; i < 6; ++ i ) {
    for( int j = 0; j < nrow; ++ j ) {
      if( yln_i[ i * nrow + j ] == 1 ) {
        yln[ i * nrow + j ] = yln_[ i * nrow + j ] - mu;
      } else {
        yln[ i * nrow + j ] = 0;
      }
    }
  }
  return;
}