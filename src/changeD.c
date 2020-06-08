#include <R.h>
#include <Rinternals.h>

SEXP changeD( SEXP old, SEXP a, SEXP b ) {
  int *lbl = INTEGER( old );
  int a_ = INTEGER( a )[ 0 ];
  int b_ = INTEGER( b )[ 0 ];
  
  int len = length( old );
  
  SEXP ans = PROTECT( allocVector( INTSXP, len ) );
  SEXP dim = getAttrib( old, R_DimSymbol );
  setAttrib( ans, R_DimSymbol, dim );
  int *ptr_ans = INTEGER( ans );
  
  for( int i = 0; i < len; ++ i ) {
    if( lbl[ i ] == a_ || lbl[ i ] == b_ ) {
      ptr_ans[ i ] = NA_INTEGER;
    } else {
      ptr_ans[ i ] = lbl[ i ];
    }
  }
  UNPROTECT( 1 );
  return ans;
}
