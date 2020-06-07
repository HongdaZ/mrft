#include <R.h>
#include <Rinternals.h>

#include <queue>
#include <stack>
#include <vector>
#include <set>

#include "search.h"
#include "findRegion.h"
#include "helper.h"
#include "initRegion.h"
#include "scTrn.h"

using std::stack;
using std::queue;
using std::vector;
using std::set;

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
    set<int> tumor_labels;
    set<int> outl_labels;
    map<int, vector<double>> health_parm;
    map<int, vector<double>> tumor_parm;
    map<int, vector<double>> outl_parm;

    
    map<int, set<int>> tumor_regions;
    initRegion( ptr_seg, ptr_nidx, len,
                tumor_regions, tumor_labels );
    
    list< map<int, int>> regions;
 
    // // Debug
    // ptr_seg[ 2 * ( 1032015 - 1 ) ] = 1;
    // ///////////////////////////////////////////////////
    int flag = scTrn( regions, tumor_labels, tumor_regions, ptr_seg,
                      ptr_nidx, 1032015 );

    // // Debug
    // int n[ 4 ] = { -32, -31, -30, -29 };
    // for( int i = 0; i < 4; ++ i ) {
    //   set<int> tmp = tumor_regions[ n[ i ] ];
    //   for( set<int>::iterator it = tmp.begin(); it != tmp.end(); ++ it ) {
    //     Rprintf( "regions: %d;", *it );
    //   }
    //   Rprintf( "\n");
    // }
    // ///////////////////////////////////////////////////
    // // Debug
    // Rprintf( "combine or split: %d \n", flag );
    // for( list<map<int, int>>::iterator it = regions.begin();
    //      it != regions.end(); ++ it ) {
    //   map<int, int> tmp_region = *it;
    //   for( map<int, int>::iterator itr = tmp_region.begin();
    //        itr != tmp_region.end(); ++ itr ) {
    //     Rprintf( "index = %d, label = %d;", itr->first, itr->second );
    //   }
    //   Rprintf("\n");
    // }
    // ///////////////////////////////////////////////////////////////
    return seg;
  }
} // extern "C"
