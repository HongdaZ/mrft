#include <R.h>
#include <Rinternals.h>

#include <list>
#include <queue>
#include <stack>

#include "search.h"
#include "findRegion.h"
#include "helper.h"

using std::stack;
using std::queue;
using std::list;

extern "C" {
  SEXP estParm( SEXP model, SEXP delta, SEXP gamma, 
                SEXP alpha, SEXP beta, SEXP lambda, 
                SEXP a, SEXP b, SEXP m, SEXP nu ) {
    SEXP info = getListElement( model, "info" );
    SEXP seg = getListElement( model, "seg" );
    
    int *ptr_label = INTEGER( seg );
    
    SEXP idx = getListElement( info, "idx" );
    SEXP nidx = getListElement( info, "nidx" );
    SEXP intst = getListElement( info, "intst" );
    SEXP nintst = getListElement( info, "nintst" );
    
    int *ptr_nidx = INTEGER( nidx );
    
    int start = 433655;
    list<int> region = findRegion( ptr_label, ptr_nidx, start );
    return nintst;
  }
} // extern "C"
