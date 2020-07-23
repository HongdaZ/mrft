#include <R.h>
#include <Rinternals.h>

extern "C" SEXP changeC( SEXP old, SEXP image, SEXP q );

SEXP changeC( SEXP old, SEXP image, SEXP q ) {
  int *lbl = INTEGER( old );
  double *img = REAL( image );
  double q_ = REAL( q )[ 0 ];

  int len = length( old );
  
  SEXP ans = PROTECT( allocVector( INTSXP, len ) );
  SEXP dim = getAttrib( old, R_DimSymbol );
  setAttrib( ans, R_DimSymbol, dim );
  int *ptr_ans = INTEGER( ans );
  
  for( int i = 0; i < len; ++ i ) {
    if( img[ i ] < q_ ) {
      ptr_ans[ i ] = NA_INTEGER;
    } else {
      ptr_ans[ i ] = lbl[ i ];
    }
  }
  UNPROTECT( 1 );
  return ans;
}
