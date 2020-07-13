#include "debug.h"

void debug() {
  vector<double> theta;
  for( int i = 0; i < 6; ++ i ) {
    theta.push_back( 0 );
  }
  int i = 0;
  for( vector<double>::iterator it = theta.begin(); 
       it != theta.end(); ++ it, ++ i ) {
    Rprintf( "theta[%d] = %f;\t", i, *it );
  }
  Rprintf( "\n" );
  vector<double> vec( 6, 0 );
  i = 0;
  for( vector<double>::iterator it = vec.begin(); 
       it != vec.end(); ++ it, ++ i ) {
    Rprintf( "vec[%d] = %f;\t", i, *it );
  }
  Rprintf( "\n" );
  return;
  
}