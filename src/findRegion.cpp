#include <R.h>
#include <Rinternals.h>
#include "search.h"
extern "C" {
  SEXP findRegion( SEXP label, SEXP nidx, SEXP start ) {
    int *ptr_label = INTEGER( label  );
    int *ptr_nidx = INTEGER(  nidx );
    int *ptr_star = INTEGER( start );
    
  }
} // extern "C"
