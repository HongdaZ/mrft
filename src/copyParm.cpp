#include "copyParm.h"

// copy estimated parameters to SEXP
void copyParm( map<int, vector<double>> &health_parm,
               map<int, vector<double>> &tumor_parm,
               map<int, vector<double>> &outl_parm,
               double *ptr_res_parm, int nrow ) {
  int j = 0;
  for( map<int, vector<double>>::iterator it = tumor_parm.begin();
       it != tumor_parm.end(); ++ it ) {
    ptr_res_parm[ nrow * j ] = it->first;
    for( int i = 0; i < 8; ++ i ) {
      ptr_res_parm[ nrow * j + i + 1 ] = it->second[ i ];
    }
    ++ j;
  }
  for( map<int, vector<double>>::iterator it = health_parm.begin();
       it != health_parm.end(); ++ it ) {
    ptr_res_parm[ nrow * j ] = it->first;
    for( int i = 0; i < 8; ++ i ) {
      ptr_res_parm[ nrow * j + i + 1 ] = it->second[ i ];
    }
    ++ j;
  }
  for( map<int, vector<double>>::iterator it = outl_parm.begin();
       it != outl_parm.end(); ++ it ) {
    ptr_res_parm[ nrow * j ] = it->first;
    for( int i = 0; i < 2; ++ i ) {
      ptr_res_parm[ nrow * j + i + 1 ] = it->second[ i ];
    }
    for( int i = 2; i < 8; ++ i ) {
      ptr_res_parm[ nrow * j + i + 1 ] = 0;
    }
    ++ j;
  }
  return;
}

void copyParmHealth( map<int, vector<double>> &health_parm,
                     double *ptr_res_parm, int nrow ) {
  int j = 0;
  for( map<int, vector<double>>::iterator it = health_parm.begin();
       it != health_parm.end(); ++ it ) {
    ptr_res_parm[ nrow * j ] = it->first;
    for( int i = 0; i < 8; ++ i ) {
      ptr_res_parm[ nrow * j + i + 1 ] = it->second[ i ];
    }
    ++ j;
  }
  return;
}
