#include "copyParm.h"

// copy estimated parameters to SEXP
void copyParm( vector<double> &health_parm,
               vector<double> &tumor_parm,
               vector<double> &outl_parm,
               double *ptr_res_parm, const int &nrow, 
               const vector<int> &tumor_labels, 
               const vector<int> &outl_labels, const int &len ) {
  int j = 0;
  int curr_label;
  for( int i = 0; i < len; ++ i ) {
    if( tumor_labels[ i ] != 0 ) {
      curr_label = - i - 4;
      ptr_res_parm[ nrow * j ] = curr_label;
      for( int k = 0; k < 8; ++ k ) {
        ptr_res_parm[ nrow * j + k + 1 ] = tumor_parm[ 8 * i + k ];
      }
      ++ j;
    } 
  }
  for( int i = 0; i < 3; ++ i ) {
    curr_label = - i - 1;
    ptr_res_parm[ nrow * j ] = curr_label;
    for( int k = 0; k < 8; ++ k ) {
      ptr_res_parm[ nrow * j + k + 1 ] = health_parm[ 8 * i + k ];
    }
    ++ j;
  }
  for( int i = 0; i < len; ++ i ) {
    if( outl_labels[ i ] != 0 ) {
      curr_label = i + 1;
      ptr_res_parm[ nrow * j ] = curr_label;
      for( int k = 0; k < 2; ++ k ) {
        ptr_res_parm[ nrow * j + k + 1 ] = outl_parm[ 2 * i + k ];
      }
      for( int k = 2; k < 8; ++ k ) {
        ptr_res_parm[ nrow * j + k + 1 ] = 0;
      }
      ++ j;
    }
  }
  return;
}
// health_parm: -1, -2, and -3
void copyParmHealth( vector<double> &health_parm,
                     double *ptr_res_parm, const int &nrow ) {
  int j = 0;
  int curr_label;
  for( int i = 0; i < 3; ++ i ) {
    curr_label = - i - 1;
    ptr_res_parm[ nrow * j ] = curr_label;
    for( int k = 0; k < 8; ++ k ) {
      ptr_res_parm[ nrow * j + k + 1 ] = health_parm[ 8 * ( 2 - i ) + k ];
    }
    ++ j;
  }
  return;
}
