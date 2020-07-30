#include <R.h>
#include <Rinternals.h>

#include <vector>
#include <set>

#include "helper.h"
#include "initParm.h"
#include "updateBeta.h"
#include "debug.h"
#include "cmpE.h"
#include "copyParm.h"
#include "restoreImg.h"

using std::vector;
using std::set;

// estimate parameters of t1ce and t2 images without tumor
extern "C" SEXP est3( SEXP model, SEXP delta, SEXP gamma, 
                        SEXP alpha, SEXP beta, SEXP lambda2, 
                        SEXP m, SEXP nu2, SEXP maxit );
SEXP est3( SEXP model, SEXP delta, SEXP gamma, 
           SEXP alpha, SEXP beta, SEXP lambda2, 
           SEXP m, SEXP nu2, SEXP maxit ) {
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
  const double *ptr_m = REAL( m );
  const double *ptr_nu2 = REAL( nu2 );
  const int *ptr_maxit = INTEGER( maxit );
  
  int len = length( idx );
  // a copy of seg and beta
  SEXP res_seg = PROTECT( allocMatrix( INTSXP, 2, len ) );
  SEXP res_beta = PROTECT( allocVector( REALSXP, 3 ) );
  int *ptr_res_seg = INTEGER( res_seg );
  for( int i = 0; i < 2 * len; ++ i ) {
    ptr_res_seg[ i ] = old_seg[ i ];
  }
  double *ptr_res_beta = REAL( res_beta );
  for( int i = 0; i < 3; ++ i ) {
    ptr_res_beta[ i ] = old_beta[ i ];
  }
  map<int, vector<double>> health_parm;
  // Initialize parameters for healthy regions
  initParmHealth3( health_parm, ptr_res_seg, ptr_m, ptr_nu2, ptr_intst, 
                   ptr_lambda2, ptr_nidx, ptr_nintst, ptr_alpha, 
                   ptr_res_beta, len, 20 );
  updateBeta3( ptr_res_beta, ptr_alpha, health_parm );
  vector<double> theta( 6, 0 );
  
  for( int i = 0; i < *ptr_maxit; ++ i ) {
    for( int j = 1; j <= len; ++ j ) {
      cmpE3( j, health_parm, ptr_res_seg,ptr_nidx, ptr_intst, ptr_nintst, 
             ptr_delta, ptr_gamma, theta );
     
    }
    initParmHealth3( health_parm, ptr_res_seg, ptr_m, ptr_nu2, ptr_intst, 
                     ptr_lambda2, ptr_nidx, ptr_nintst, ptr_alpha, 
                     ptr_res_beta, len, 20 );
  }
  // segment zero blocks
  for( int j = 1; j <= len; j ++ ) {
    if( ptr_res_seg[ 2 * ( j - 1 ) ] == 0 ) {
      cmpE3( j, health_parm, ptr_res_seg,ptr_nidx, ptr_intst, ptr_nintst, 
             ptr_delta, ptr_gamma, theta );
    }
  }
  
  int nrow = 1 + 2 + 6;
  int ncol = health_parm.size();
  
  SEXP res_parm = PROTECT( allocMatrix( REALSXP, nrow, ncol ) );
  SEXP res_image = PROTECT( alloc3DArray( INTSXP, 240, 240, 155 ) );
  double *ptr_res_parm = REAL( res_parm );
  int *ptr_res_image = INTEGER( res_image );
  copyParmHealth( health_parm, ptr_res_parm, nrow );
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