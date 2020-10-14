#include <R.h>
#include <Rinternals.h>

#include <vector>
#include <list>

#include "helper.h"
#include "cnctRegion.h"
#include "excldVoxel.h"
#include "extRegion.h"
#include "pad2zero.h"

using std::vector;
using std::list;

// Postprocess the results
extern "C" SEXP postProcess( SEXP post_data );
SEXP postProcess( SEXP post_data ) {
  SEXP t1ce = getListElement( post_data, "t1ce_seg" );
  SEXP flair = getListElement( post_data, "flair_seg" );
  SEXP t2 = getListElement( post_data, "t2_seg" );
  SEXP idx = getListElement( post_data, "idx" );
  SEXP nidx = getListElement( post_data, "nidx" );
  SEXP aidx = getListElement( post_data, "aidx" );
  
  int *ptr_t1ce = INTEGER( t1ce );
  int *ptr_flair = INTEGER( flair );
  int *ptr_t2 = INTEGER( t2 );
  const int *ptr_idx = INTEGER( idx );
  const int *ptr_nidx = INTEGER( nidx );
  const int *ptr_aidx = INTEGER( aidx );
  
  const int len = length( idx );
  SEXP res = PROTECT( allocMatrix( INTSXP, 2, len ) );
  // Store the results for each tissue type
  int *ptr_res = INTEGER( res );
  int *ptr_hemorrhage = new int[ 2 * len ]();
  int *ptr_necrosis = new int[ 2 * len ]();
  int *ptr_nonenh = new int[ 2 * len ]();
  int *ptr_enh = new int[ 2 * len ]();
  int *ptr_edema = new int[ 2 * len ]();
  
  // Store the result of findRegion
  vector<int> region;
  region.reserve( len );
  
  // 10-1: Find hemorrhage
  for( int i = 0; i < len; ++ i ) {
    if( ptr_flair[ 2 * i ] == 1 && 
        ptr_t2[ 2 * i ] == 1 &&
        ptr_t1ce[ 2 * i ] != 4 ) {
      ptr_hemorrhage[ 2 * i ] = 5;
    }
  }
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_hemorrhage, ptr_flair, 
                    1, region ) ) {
      excldVoxel( region, ptr_t2, 4 );
      excldVoxel( region, ptr_t1ce, 4 );
      extRegion( region, ptr_hemorrhage, 5, .5 );
    } 
  }
  // Recover the padding to zero
  pad2zero( ptr_hemorrhage, len );
  
  // 10-2: Find necrosis
  for( int i = 0; i < len; ++ i ) {
    if( ( ptr_t1ce[ 2 * i ] == 1 || ptr_flair[ 2 * i ] == 1 ) &&
        ptr_flair[ 2 * i ] != 4 ) {
      ptr_necrosis[ 2 * i ] = 6;
    }
  }
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_necrosis, ptr_t2, 
                    4, region ) ) {
      excldVoxel( region, ptr_t1ce, 4 );
      excldVoxel( region, ptr_flair, 4 );
      excldVoxel( region, ptr_hemorrhage, 5 );
      extRegion( region, ptr_necrosis, 6, 0 );
    }
  }
  // Recover the padding to zero
  pad2zero( ptr_necrosis, len );
  
  // 10-3: Find enhancing tumor core
  for( int i = 0; i < len; ++ i ) {
    if( ( ptr_t1ce[ 2 * i ] == 4 && ptr_necrosis[ 2 * i ] != 6 ) &&
        ( ptr_t2[ 2 * i ] == 4 || ptr_flair[ 2 * i ] == 4 ) ) {
      ptr_enh[ 2 * i ] = 4;
    }
  }
  
  // 10-4: Find edema
  for( int i = 0; i < len; ++ i ) {
    if( ptr_flair[ 2 * i ] == 4 || ptr_t2[ 2 * i ] == 4 ||
        ptr_enh[ 2 * i ] == 4 && ptr_necrosis[ 2 * i ] != 6 &&
        ptr_hemorrhage[ 2 * i ] != 5 ) {
      ptr_edema[ 2 * i ] = 2;
    }
  }
  
  
  delete [] ptr_hemorrhage;
  delete [] ptr_necrosis;
  delete [] ptr_nonenh;
  delete [] ptr_enh;
  delete [] ptr_edema;
  UNPROTECT( 1 );
  return res;
}