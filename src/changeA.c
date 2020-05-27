#include <R.h>
#include <Rinternals.h>

SEXP changeA( SEXP img, SEXP label ) {
  
  double *image = REAL( img );
  int *lbl = INTEGER( label );
  int len = length( label );
  SEXP ans = PROTECT( allocVector( REALSXP, len ) );
  double *ptr_ans = REAL( ans );
  for( int i = 0; i < len; ++ i ) {
    if( lbl[ i ] == 0 ) {
      ptr_ans[ i ] = image[ i ];
    } else {
      ptr_ans[ i ] = R_NaN;
    }
  }
  UNPROTECT( 1 );
  return ans;
}
