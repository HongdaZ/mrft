#include <R.h>
#include <Rinternals.h>

SEXP changeB( SEXP label, SEXP image1, SEXP q1, SEXP image2, SEXP q2, 
              SEXP k  ) {
  int *lbl = INTEGER( label );
  double *img1 = REAL( image1 );
  double *img2 = REAL( image2 );
  double q_1 = REAL( q1 )[ 0 ];
  double q_2 = REAL( q2 )[ 0 ];
  int l = INTEGER( k )[ 0 ];
  int len = length( label );

  for( int i = 0; i < len; ++ i ) {
    if( img1[ i ] < q_1 && img2[ i ] > q_2 ) {
      lbl[ i ] = l;
    }
  }

  return R_NilValue;
}
