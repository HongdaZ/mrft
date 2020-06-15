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
#include "updateMu.h"
#include "updateTS.h"
#include "updateSigma.h"
#include "energyY.h"

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
    // // For testing scTrn
    // int flag = scTrn( regions, tumor_labels, tumor_regions, ptr_seg,
    //                   ptr_nidx, 1032015 );
    // Rprintf( "combine or split: %d \n", flag );
    // ///////////////////////////////////////////////////////////////
    // // Debug updateMu
    // for( int j = 1; j < 4; ++ j ) {
    //   set<int> region_h;
    //   for( int i = 0; i < len; ++ i ) {
    //     if( ptr_seg[ 2 * i ] == - j ) {
    //       region_h.insert( i + 1 ); // region starts from 1
    //     }
    //   }
    //   double sigma2 = 2;
    //   double m = 3;
    //   double nu2 = 2;
    //   vector<double> theta;
    //   for( int k = 0; k < 6; ++ k ) {
    //     theta.push_back( .01 + k * .01 );
    //   }
    //   double mu = updateMu( region_h, sigma2, m, nu2, theta, ptr_intst );
    //   Rprintf( " region %d, mu = %f", j, mu );
    // }
    //////////////////////////////////////////////////////////////////
    // // Debug updateMu for tumor
    // map<int, int> tumor_34;
    // for( int i = 0; i < len; ++ i ) {
    //   if( ptr_seg[ 2 * i ] == -34 ) {
    //     tumor_34[ i + 1 ] = -34;
    //   }
    // }
    // vector<double> theta;
    // for( int k = 0; k < 6; ++ k ) {
    //   theta.push_back( .01 + k * .01 );
    // }
    // double mu = updateMu( tumor_34, 2, 1, 2, 5, 6, theta, ptr_intst );
    // Rprintf( " Tumor region -34, mu = %f", mu );
    // // Debug updateTS
    // set<int> region_h;
    // for( int i = 0; i < len; ++ i ) {
    //   if( ptr_seg[ 2 * i ] == - 23 ) {
    //     region_h.insert( i + 1 ); // region starts from 1
    //   }
    // }
    // double sigma2 = 2;
    // vector<double> theta;
    // updateTS( region_h, -23, 1, sigma2, 3, ptr_seg, ptr_nidx, ptr_intst, 
    //           ptr_nintst, theta, 1, .00001 );
    // // updateSigma( 20, 1, ptr_intst, 2, .5, sigma2 );
    // // Rprintf( "sigma2 = %f\n", sigma2 );
    // // for( int k = 0; k < 6; ++ k ) {
    // //   Rprintf( "theta = %f\t", theta[ k ] );
    // // }
    // // Rprintf( "\n" );
    // //////////////////////////////////////////////////////////////////////////
    // double energy = energyY( region_h, 1, .5, sigma2, 3, ptr_seg, ptr_nidx,
    //                          ptr_intst, ptr_nintst, theta, 
    //                          1, .00001, ptr_a[ 0 ], ptr_b[ 0 ] );
    // // Rprintf( "energyY = %f\n", energy );
    vector<double> theta;
    double vtheta[ 6 ] = { 0.1, 0.3, 0.4, 0.2, 0.5, 6.0 };
    for( int i = 0; i < 6; ++ i ) {
      theta.push_back( vtheta[ i ] ) ;
    }
    double energy = energyY( -1, 1208, 1, 0.13, ptr_seg, 
                             ptr_nidx, ptr_intst, ptr_nintst, theta );
    Rprintf( "energyY = %f\n", energy);
    return seg;
  }
} // extern "C"
