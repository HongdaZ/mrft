#include "debug.h"

// void debug() {
//   vector<double> theta;
//   for( int i = 0; i < 6; ++ i ) {
//     theta.push_back( 0 );
//   }
//   int i = 0;
//   for( vector<double>::iterator it = theta.begin(); 
//        it != theta.end(); ++ it, ++ i ) {
//     Rprintf( "theta[%d] = %f;\t", i, *it );
//   }
//   Rprintf( "\n" );
//   vector<double> vec( 6, 0 );
//   i = 0;
//   for( vector<double>::iterator it = vec.begin(); 
//        it != vec.end(); ++ it, ++ i ) {
//     Rprintf( "vec[%d] = %f;\t", i, *it );
//   }
//   Rprintf( "\n" );
//   return;
//   
// }
void debug() {
  int curr_label = 1;
  double mu = 2;
  double sigma2 = 3;
  vector<double> theta = { 9, 8, 7, 6, 5, 4 };
  vector<double> tmp_parm; 
  tmp_parm.push_back( curr_label );
  tmp_parm.push_back( mu );
  tmp_parm.push_back( sigma2 );
  tmp_parm.insert( tmp_parm.end(), theta.begin(), theta.end() );
  
  int i = 0;
  for( vector<double>::iterator it = tmp_parm.begin();
       it != tmp_parm.end(); ++ it, ++ i ) {
    Rprintf( "tmp_parm[%d] = %f;\t", i, *it );
  }
  Rprintf( "\n" );
  
  vector<double> tap_parm( 9, 0 ); 
  tap_parm[ 0 ] = curr_label;
  tap_parm[ 1 ] = mu;
  tap_parm[ 2 ] = sigma2;
  for( i = 0; i < 6; ++ i ) {
    tap_parm[ i + 3 ] = theta[ i ];
  }
  i = 0;
  for( vector<double>::iterator it = tap_parm.begin();
       it != tap_parm.end(); ++ it, ++ i ) {
    Rprintf( "tap_parm[%d] = %f;\t", i, *it );
  }
  Rprintf( "\n" );
  return;
}