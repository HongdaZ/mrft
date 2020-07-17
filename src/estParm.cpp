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
#include "scPred.h"
#include "updateMu.h"
#include "updateTS.h"
#include "updateSigma.h"
#include "energyY.h"
#include "energyX.h"
#include "nbrLabel.h"
#include "updateParm.h"
#include "initParm.h"
#include "cmpET.h"
#include "cmpEP.h"
#include "updateBeta.h"
#include "skip.h"
#include "debug.h"
#include "updateHealth.h"
#include "updateTumor.h"

using std::stack;
using std::queue;
using std::vector;
using std::set;

extern "C" {
  SEXP estParm( SEXP model, SEXP delta, SEXP gamma, 
                SEXP alpha, SEXP beta, SEXP lambda2, 
                SEXP a, SEXP b, SEXP m, SEXP nu2 ) {
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
    double *ptr_beta = REAL( beta );
    const double *ptr_lambda2 = REAL( lambda2 );
    const double *ptr_a = REAL( a );
    const double *ptr_b = REAL( b );
    const double *ptr_m = REAL( m );
    const double *ptr_nu2 = REAL( nu2 );
    
    int len = length( idx );
    list<int> tumor_labels;
    list<int> outl_labels;
    map<int, vector<double>> health_parm;
    map<int, vector<double>> tumor_parm;
    map<int, vector<double>> outl_parm;

    
    map<int, list<int>> tumor_regions;
    initRegion( ptr_seg, ptr_nidx, len,
                tumor_regions, tumor_labels );
    Rprintf( "initRegion finished!\n" );
    // Rprintf( "tumor_regions.size = %d; tumor_labels.size() = %d\n", 
    //          tumor_regions.size(), tumor_labels.size() );
    // list<int>::iterator it_label = tumor_labels.begin();
    // map<int, list<int>>::iterator it_region = tumor_regions.begin();
    // for( ;  it_label != tumor_labels.end(); ++ it_label, ++ it_region ) {
    //   Rprintf( "label = %d; size of region = %d\n", it_region->first, 
    //            it_region->second.size() );
    // }
    // Rprintf( "label = %d; idx = %d", -15, tumor_regions[ - 15 ].front() );
    int curr_label = - 4;
    int maxit = 20;
    list<int> &region = tumor_regions[ - 4 ];
    double mu = -1, sigma2 = 1; // sigma2 has to be non-zero;
    vector<double> theta( 6, 0 );
    updateParm( mu, theta, sigma2, region, ptr_m[ 3 ], ptr_m[ 2 ],
                ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst, curr_label, 
                ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, ptr_nintst,
                ptr_alpha[ 3 ], ptr_beta[ 3 ], maxit );
    vector<double> t_parm( 8, 0 );
    t_parm[ 0 ] = mu;
    t_parm[ 1 ] = sigma2;
    for( int i = 0; i < 6; ++ i ) {
      t_parm[ i + 2 ] = theta[ i ];
    }
    // tumor_parm[ curr_label ] = t_parm;
    // Rprintf( "label = %d; mu = %f; sigma2 = %f; theta = ", curr_label, mu,
    //          sigma2 );
    // for( int i = 0; i < 6; ++ i ) {
    //   Rprintf( "%f, ", theta[ i ] );
    // }
    // Rprintf( "\n" );
    // initParm( health_parm, tumor_parm, ptr_seg, ptr_m, ptr_nu2, ptr_intst,
    //           ptr_lambda2, ptr_nidx, ptr_nintst, ptr_alpha, ptr_beta,
    //           tumor_regions, ptr_a, ptr_b, len, 20 );
    // Rprintf( "initParm finished!\n" );
    // updateBeta( ptr_beta, ptr_alpha, health_parm, tumor_regions,
    //             tumor_parm );
    // Rprintf( "updateBeta finished!\n" );
    
    // int maxit = 20;
    // bool skip_curr;
    // vector<int> search( len, 0 );
    // double threshold = ptr_m[ 2 ];
    // Rprintf( "Segmentation started!\n" );
    // for( int i = 0; i < maxit; ++ i ) {
    //   for( int j = 1; j <= len; ++ j ) {
    //     
    //     skip_curr = skip( j, ptr_seg, ptr_nidx, ptr_intst, threshold, 
    //                       search[ j ], 10 );
    //     if( skip_curr ) {
    //       continue;
    //     } else {
    //       ++ search[ j ];
    //       list<list<int>> regions;
    //       list<int> labels;
    //       int sc = scTrn( labels, regions, tumor_labels, tumor_regions, 
    //                       ptr_seg, ptr_nidx, j );
    //       cmpET( j, sc, labels, regions, tumor_regions, tumor_labels, 
    //              outl_labels, health_parm, tumor_parm, outl_parm, ptr_seg,
    //              ptr_nidx, ptr_intst, ptr_nintst, ptr_delta, ptr_gamma,
    //              ptr_alpha, ptr_beta, ptr_lambda2, ptr_a, ptr_b, ptr_m,
    //              ptr_nu2 );
    //     }
    //     Rprintf( "%d\t; curr_label = %d\n", j, ptr_seg[ 2 * ( j - 1 ) ] );
    //   }
    //   // update parm for healthy regions
    //   updateHealth( health_parm, ptr_seg, len, ptr_m, ptr_nu2, ptr_intst,
    //                 ptr_lambda2, ptr_nidx, ptr_nintst, ptr_alpha, ptr_beta,
    //                 maxit );
    //   // update parm for tumor_regions
    //   updateTumor( tumor_labels, tumor_parm, tumor_regions, ptr_a, ptr_b,
    //                ptr_seg, ptr_m, ptr_intst, ptr_lambda2, ptr_nidx,
    //                ptr_nintst, ptr_alpha, ptr_beta, maxit );
    // }
    return seg;
  }
} // extern "C"
