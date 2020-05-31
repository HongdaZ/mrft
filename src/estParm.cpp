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
    
    int *ptr_idx = INTEGER( idx );
    int *ptr_nidx = INTEGER( nidx );
    double *ptr_intst = REAL( intst );
    double *ptr_nintst = REAL( nintst );
    
    int len = length( idx );
    list<int> tumor_labels;
    list<int> out_labels;
    
    map<int, list<int>> tumor_regions = initRegion( ptr_seg, ptr_nidx, len,
                                                    tumor_labels );
    return seg;
  }
} // extern "C"
