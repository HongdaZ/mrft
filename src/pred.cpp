#include <R.h>
#include <Rinternals.h>

#include <vector>
#include <set>
#include <chrono> 

#include "helper.h"
#include "initRegion.h"
#include "scPred.h"
#include "initParm.h"
#include "cmpEP.h"
#include "updateBeta.h"
#include "skip.h"
#include "debug.h"
#include "copyParm.h"
#include "restoreImg.h"

using std::vector;
using std::set;
using namespace std::chrono; 

// predict the labels for T1ce or T2 images
extern "C" SEXP pred4( SEXP model, SEXP delta, SEXP gamma, 
                      SEXP alpha, SEXP beta, SEXP lambda2, 
                      SEXP a, SEXP b, SEXP m, SEXP nu2, SEXP maxit );

SEXP pred4( SEXP model, SEXP delta, SEXP gamma, 
            SEXP alpha, SEXP beta, SEXP lambda2, 
            SEXP a, SEXP b, SEXP m, SEXP nu2, SEXP maxit ) {
  SEXP info = getListElement( model, "info" );
  SEXP seg = getListElement( model, "seg" );
  
  const int *old_seg = INTEGER( seg );
  
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
  const double *old_beta = REAL( beta );
  const double *ptr_lambda2 = REAL( lambda2 );
  const double *ptr_a = REAL( a );
  const double *ptr_b = REAL( b );
  const double *ptr_m = REAL( m );
  const double *ptr_nu2 = REAL( nu2 );
  const int *ptr_maxit = INTEGER( maxit );
  
  int len = length( idx );
  // a copy of seg and beta
  SEXP res_seg = PROTECT( allocMatrix( INTSXP, 2, len ) );
  SEXP res_beta = PROTECT( allocVector( REALSXP, 4 ) );
  int *ptr_res_seg = INTEGER( res_seg );
  for( int i = 0; i < 2 * len; ++ i ) {
    ptr_res_seg[ i ] = old_seg[ i ];
  }
  double *ptr_res_beta = REAL( res_beta );
  for( int i = 0; i < 4; ++ i ) {
    ptr_res_beta[ i ] = old_beta[ i ];
  }
  vector<int> tumor_labels( len, 0 );
  vector<int> outl_labels( len, 0 );
  int n_tumor = 0;
  int n_outl = 0;
  // 0, 1, 2 = -1, -2, -3
  vector<double> health_parm( 8 * 3, 0 );
  // 0, 1, 2, ... = -4, -5, -6, ...
  vector<double> tumor_parm( 8 * len, 0 );
  // 0, 1, 2, ... = 1, 2, 3, ...
  vector<double> outl_parm( 2 * len, 0 );
  list<list<int>> tumor_regions;
  
  // whole and subregions in sc functions
  vector<int> regions_whole;
  regions_whole.reserve( len );
  vector<int> regions_sub;
  regions_sub.reserve( len );
  
  // frontier in search
  vector<int> front;
  front.reserve( len );
  // store the result of findRegion
  vector<int> region;
  region.reserve( len );
  
  // number labels
  vector<int> labels;
  labels.reserve( 6 );
  
  bool skip_curr;
  vector<int> search( len, 0 );
  double lower = ptr_m[ 2 ];
  double upper = ptr_m[ 3 ];
  vector<double> outlier_parm( 3, 0 );
  vector<double> theta( 6, 0 );
  vector<double> tmp_parm( 9, 0 );
  vector<double> out_theta( 6, 0 );
  vector<double> new_out_parm( 2, 0 );
  vector<double> whole_parm( 8, 0 );
  
  initRegion( region, front, tumor_regions,
              ptr_res_seg, ptr_nidx, len,
              tumor_labels, n_tumor );
  Rprintf( "n_tumor = %d\n", n_tumor );
  initParm( region, theta, true, health_parm, tumor_parm, ptr_res_seg,
            ptr_m, ptr_nu2,
            ptr_intst, ptr_lambda2, ptr_nidx, ptr_nintst, ptr_alpha,
            ptr_res_beta, tumor_regions,
            n_tumor, ptr_a, ptr_b, len, 20 );
  updateBeta( ptr_res_beta, ptr_alpha, health_parm, n_tumor,
              tumor_parm );
  int old_label = 0;
  int new_label = 0;
  int curr_idx = 0;
  auto start = high_resolution_clock::now(); 
  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start); 
  int spend = 0; 
  int max_idx;
  for( int i = 0; i < *ptr_maxit; ++ i ) {
    for( int j = 1; j <= len; ++ j ) {
      curr_idx = j;
      // skip the voxels whose label remain the same in 5 consecutive
      // updates
      // skip_curr = skip( j, ptr_res_seg, ptr_nidx, ptr_intst, lower, upper,
      //                   search[ j - 1 ], 3 );

      old_label = ptr_res_seg[ 2 * ( curr_idx - 1 ) ];
      // Get starting timepoint
      start = high_resolution_clock::now();
      int sc = scPred( labels, regions, front, region, tumor_labels,
                       n_voxel, ptr_res_seg, ptr_nidx, len, curr_idx );
      // Get ending timepoint
      stop = high_resolution_clock::now();
      // Get duration. Substart timepoints to
      // get durarion. To cast it to proper unit
      // use duration cast method
      duration = duration_cast<microseconds>(stop - start);
      if( spend < duration.count() ) {
        spend =  duration.count();
        max_idx = curr_idx;
      }

      // cmpEP( region, curr_idx, sc, labels, regions, tumor_labels,
      //        outl_labels, health_parm, tumor_parm, outl_parm, ptr_res_seg,
      //        ptr_nidx, ptr_intst, ptr_nintst, ptr_delta, ptr_gamma,
      //        ptr_alpha, ptr_res_beta, ptr_lambda2, ptr_a, ptr_b, ptr_m,
      //        ptr_nu2, outlier_parm, theta, tmp_parm, out_theta,
      //        new_out_parm, whole_parm, n_tumor, n_outl, len );
      new_label = ptr_res_seg[ 2 * ( curr_idx - 1 ) ];
      // Rprintf( "j = %d,\t curr_idx = %d,\t sc = %d,\t old = %d,\t new = %d,\t time = %d\n",
      //          j, curr_idx, sc, old_label, new_label, duration.count() );
    }
    Rprintf( "update parm for healthy and tumorous\n" );
    //update parm for healthy and tumorous regions
    initParm( region, theta, false, health_parm, tumor_parm, ptr_res_seg, ptr_m, ptr_nu2,
              ptr_intst, ptr_lambda2, ptr_nidx, ptr_nintst, ptr_alpha,
              ptr_res_beta, n_voxel, n_tumor, ptr_a, ptr_b, len, 20 );
  }
  // segment zero blocks
  for( int j = 1; j <= len; j ++ ) {
    if( ptr_res_seg[ 2 * ( j - 1 ) ] == 0 ) {
      int curr_idx = j;
      int sc = scPred( labels, regions, front, region, tumor_labels,
                       n_voxel, ptr_res_seg, ptr_nidx, len, curr_idx );
      cmpEP( region, curr_idx, sc, labels, regions, tumor_labels,
             outl_labels, health_parm, tumor_parm, outl_parm, ptr_res_seg,
             ptr_nidx, ptr_intst, ptr_nintst, ptr_delta, ptr_gamma,
             ptr_alpha, ptr_res_beta, ptr_lambda2, ptr_a, ptr_b, ptr_m,
             ptr_nu2, outlier_parm, theta, tmp_parm, out_theta,
             new_out_parm, whole_parm, n_tumor, n_outl, len );
    }
  }
  int nrow = 1 + 2 + 6;
  int ncol = n_tumor + 3 + n_outl;

  SEXP res_parm = PROTECT( allocMatrix( REALSXP, nrow, ncol ) );
  SEXP res_image = PROTECT( alloc3DArray( INTSXP, 240, 240, 155 ) );
  double *ptr_res_parm = REAL( res_parm );
  int *ptr_res_image = INTEGER( res_image );
  copyParm( health_parm, tumor_parm, outl_parm, ptr_res_parm, nrow,
            tumor_labels, outl_labels, len );
  restoreImg( ptr_idx, ptr_res_seg, ptr_res_image, len );

  // results to list
  SEXP names = PROTECT( allocVector( STRSXP, 4 ) );

  SET_STRING_ELT( names, 0, mkChar( "beta" ) );
  SET_STRING_ELT( names, 1, mkChar( "seg" ) );
  SET_STRING_ELT( names, 2, mkChar( "parm" ) );
  SET_STRING_ELT( names, 3, mkChar( "image" ) );

  SEXP res = PROTECT( allocVector( VECSXP, 4 ) );
  SET_VECTOR_ELT( res, 0, res_beta );
  SET_VECTOR_ELT( res, 1, res_seg );
  SET_VECTOR_ELT( res, 2, res_parm );
  SET_VECTOR_ELT( res, 3, res_image );
  setAttrib( res, R_NamesSymbol, names );

  UNPROTECT( 6 );
  return res;
}