#include "assignParm.h"

// assign parm to tumor_parm or health_parm
void assignParm( vector<double> &x_parm, const int &curr_label,
                  const vector<double> &parm ) {
  int cidx;
  if( curr_label < -3 ) {
    cidx = - curr_label - 4;
  } else if( curr_label >= -3 && curr_label <= -1 ) {
    cidx = - curr_label - 1;
  }
  for( int j = 0; j < 8; ++ j ) {
    x_parm[ 8 * cidx + j ] = parm[ j ];
  }
}
void assignParm( vector<double> &x_parm, const int &curr_label, 
                  const double &mu, const double &sigma2,
                  const vector<double> &theta ) {
  int cidx;
  if( curr_label < -3 ) {
    cidx = - curr_label - 4;
  } else if( curr_label >= -3 && curr_label <= -1 ) {
    cidx = - curr_label - 1;
  }
  x_parm[ 8 * cidx + 0 ] = mu;
  x_parm[ 8 * cidx + 1 ] = sigma2;
  for( int j = 0; j < 6; ++ j ) {
    x_parm[ 8 * cidx + j + 2 ] = theta[ j ];
  }
}