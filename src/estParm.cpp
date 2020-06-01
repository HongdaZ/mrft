#include <R.h>
#include <Rinternals.h>

#include <list>
#include <queue>
#include <stack>

#include "search.h"
#include "findRegion.h"
#include "helper.h"
#include "initRegion.h"

using std::stack;
using std::queue;
using std::list;

extern "C" {
  SEXP estParm( SEXP model, SEXP delta, SEXP gamma, 
                SEXP alpha, SEXP beta, SEXP lambda, 
                SEXP a, SEXP b, SEXP m, SEXP nu ) {
    SEXP info = getListElement( model, "info" );
    SEXP seg = getListElement( model, "seg" );
    
    int *ptr_seg = INTEGER( seg );
    
    SEXP idx = getListElement( info, "idx" );
    SEXP nidx = getListElement( info, "nidx" );
    SEXP intst = getListElement( info, "intst" );
    SEXP nintst = getListElement( info, "nintst" );
    
    const int *ptr_idx = INTEGER( idx );
    const int *ptr_nidx = INTEGER( nidx );
    const double *ptr_intst = REAL( intst );
    const double *ptr_nintst = REAL( nintst );
    
    const double *ptr_delta = REAL( delta );
    const double *ptr_gamma = REAL( gamma );
    const double *ptr_alpha = REAL( alpha );
    const double *ptr_beta = REAL( beta );
    const double *ptr_lambda = REAL( lambda );
    const double *ptr_a = REAL( a );
    const double *ptr_b = REAL( b );
    const double *ptr_m = REAL( m );
    const double *ptr_nu = REAL( nu );
    
    int len = length( idx );
    list<int> tumor_labels;
    list<int> outl_labels;

    
    map<int, list<int>> tumor_regions = initRegion( ptr_seg, ptr_nidx, len,
                                                    tumor_labels );
    return seg;
  }
} // extern "C"
