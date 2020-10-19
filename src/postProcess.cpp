#include <R.h>
#include <Rinternals.h>

#include <vector>
#include <list>

#include "helper.h"
#include "cnctRegion.h"
#include "excldVoxel.h"
#include "extRegion.h"
#include "pad2zero.h"
#include "regions.h"
#include "region2slice.h"
#include "enclose.h"
#include "inRegion.h"
#include "excldRegion.h"

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
  SEXP r_nr = getListElement( post_data, "nr" );
  SEXP r_nc = getListElement( post_data, "nc" );
  SEXP r_ns = getListElement( post_data, "ns" );
  
  const int nr = INTEGER( r_nr )[ 0 ];
  const int nc = INTEGER( r_nc )[ 0 ];
  const int ns = INTEGER( r_ns )[ 0 ];
  
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
  int *ptr_enclose_nec = new int[ 2 * len ]();
  int *ptr_enclose_hem = new int[ 2 * len ]();
  
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
    if( ( ptr_flair[ 2 * i ] == 4 || ptr_t2[ 2 * i ] == 4 ||
        ptr_enh[ 2 * i ] == 4 ) && 
        ( ptr_necrosis[ 2 * i ] != 6 &&
        ptr_hemorrhage[ 2 * i ] != 5 ) ) {
      ptr_edema[ 2 * i ] = 2;
    }
  }
  // 10-5: Find necrosis
  // necrosis enclosed by edema
  inRegion( ptr_enclose_nec, len, ptr_edema, 2, ptr_necrosis, 6, 
            region, ptr_nidx, ptr_aidx, nr, nc, ns );
  // Remove necrosis regions separate from edema
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_enclose_nec, ptr_enclose_nec, 
                    1, region ) ) {
      excldRegion( region, ptr_nidx, ptr_enclose_nec, 
                   ptr_edema, 2 );
    }
  }
  // Recover the padding to zero
  pad2zero( ptr_enclose_nec, len );
  // Extend in ptr_necrosis
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_enclose_nec, ptr_necrosis, 
                    6, region ) ) {
      extRegion( region, ptr_enclose_nec, 1, .5 );
    }
  }
  // Recover the padding to zero
  pad2zero( ptr_enclose_nec, len );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_enclose_nec[ 2 * i ] == 1 ) {
      ptr_necrosis[ 2 * i ] = 6;
    } else {
      ptr_necrosis[ 2 * i ] = 0;
    }
  } 
  pad2zero( ptr_necrosis, len );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_necrosis[ 2 * i ] == 6 ) {
      ptr_edema[ 2 * i ] = 0;
    }
  }
  
  // 10-6: Find hemorrhage
  // hemorrhage enclosed by edema
  inRegion( ptr_enclose_hem, len, ptr_edema, 2, ptr_hemorrhage, 5,
            region, ptr_nidx, ptr_aidx, nr, nc, ns );
  // Remove hemorrhage regions separate from edema
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_enclose_hem, ptr_enclose_hem,
                    1, region ) ) {
      excldRegion( region, ptr_nidx, ptr_enclose_hem,
                   ptr_edema, 2 );
    }
  }
  // Recover the padding to zero
  pad2zero( ptr_enclose_hem, len );
  // Extend in ptr_hemorrhage
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_enclose_hem, ptr_hemorrhage,
                    5, region ) ) {
      extRegion( region, ptr_enclose_hem, 1, .5 );
    }
  }
  // Recover the padding to zero
  pad2zero( ptr_enclose_hem, len );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_enclose_hem[ 2 * i ] == 1 ) {
      ptr_hemorrhage[ 2 * i ] = 5;
    } else {
      ptr_hemorrhage[ 2 * i ] = 0;
    }
  }
  pad2zero( ptr_hemorrhage, len );
  
  // for( int i = 0; i < len; ++ i ) {
  //   ptr_res[ 2 * i ] = ptr_enclose_hem[ 2 * i ];
  // }
  
  delete [] ptr_hemorrhage;
  delete [] ptr_necrosis;
  delete [] ptr_nonenh;
  delete [] ptr_enh;
  delete [] ptr_edema;
  delete [] ptr_enclose_nec;
  delete [] ptr_enclose_hem;
  
  UNPROTECT( 1 );
  return res;
}