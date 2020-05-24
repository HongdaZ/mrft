#include <R.h>
#include <Rinternals.h>
#include "search.h"
using std::stack;
using std::queue;
using std::size_t;

extern "C" {
  SEXP findRegion( SEXP label, SEXP nidx, SEXP start ) {
    int *ptr_label = INTEGER( label  );
    int *ptr_nidx = INTEGER(  nidx );
    int *ptr_start = INTEGER( start );
    
    queue<int> front;
    stack<int> region;
    
    int l = ptr_label[ 2 * ( ptr_start[ 0 ] - 1 ) ];
    front.push( ptr_start[ 0 ] );
    region.push( ptr_start[ 0 ] );
    ptr_label[ 2 * ptr_start[ 0 ] - 1 ] = 1;
    search( ptr_label, ptr_nidx, front, l, region );
    int len = region.size();
    SEXP res = PROTECT( allocVector( INTSXP, len ) );
    int *ptr_res = INTEGER( res );
    for( int i = 0; i < len; ++ i ) {
      ptr_res[ i ] = region.top();
      ptr_label[ 2 * region.top() - 1 ] = 0;
      // Rprintf( "region = %d \n", region.top() );
      region.pop();
    }
    Rprintf( "NA_INTEGER = %d \n", NA_INTEGER );
    UNPROTECT( 1 );
    return res;
  }
} // extern "C"
