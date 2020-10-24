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
#include "restoreImg.h"
#include "tissueType.h"
#include "wrapUp.h"
#include "furtherSeg.h"

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
  SEXP res_image = PROTECT( alloc3DArray( INTSXP, nr, nc, ns ) );
  SEXP hgg = PROTECT( ScalarInteger( 0 ) );
  
  // Store the results for each tissue type
  int *ptr_seg = new int[ 2 * len ]();
  int *ptr_hgg = INTEGER( hgg );
  int *ptr_res_image = INTEGER( res_image );
  int *ptr_hemorrhage = new int[ 2 * len ]();
  int *ptr_necrosis = new int[ 2 * len ]();
  int *ptr_enh = new int[ 2 * len ]();
  int *ptr_edema = new int[ 2 * len ]();
  // Necrosis inside edema
  int *ptr_enclose_nec = new int[ 2 * len ]();
  // Hemorrhage inside edema
  int *ptr_enclose_hem = new int[ 2 * len ]();
  // Necrosis inside enh
  int *ptr_enclose_ncr = new int[ 2 * len ]();
  int *ptr_tumor = new int[ 2 * len ]();
  // Store the result of findRegion
  vector<int> region;
  region.reserve( len );
  const int min_enh = 2000;
  const int min_tumor = 20000;
  
  // 10-1: Find hemorrhage
  for( int i = 0; i < len; ++ i ) {
    if( ptr_flair[ 2 * i ] == Flair::FCSF && 
        ptr_t2[ 2 * i ] == T2::T2WM &&
        ptr_t1ce[ 2 * i ] != T1ce::T1TM ) {
      ptr_hemorrhage[ 2 * i ] = Tumor::HMG;
    }
  }
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_hemorrhage, ptr_flair, 
                    Flair::FCSF, region ) ) {
      excldVoxel( region, ptr_t2, T2::T2CSF );
      excldVoxel( region, ptr_t1ce, T1ce::T1TM );
      extRegion( region, ptr_hemorrhage, Tumor::HMG, .5 );
    } 
  }
  // Recover the padding to zero
  pad2zero( ptr_hemorrhage, len );
  
  // 10-2: Find necrosis
  for( int i = 0; i < len; ++ i ) {
    if( ( ptr_t1ce[ 2 * i ] == T1ce::T1CSF || 
        ptr_flair[ 2 * i ] == Flair::FCSF ) &&
        ptr_flair[ 2 * i ] != Flair::FTM ) {
      ptr_necrosis[ 2 * i ] = Tumor::NCR;
    }
  }
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_necrosis, ptr_t2, 
                    T2::T2CSF, region ) ) {
      excldVoxel( region, ptr_t1ce, T1ce::T1TM );
      excldVoxel( region, ptr_flair, Flair::FTM );
      excldVoxel( region, ptr_hemorrhage, Tumor::HMG );
      extRegion( region, ptr_necrosis, Tumor::NCR, 0 );
    }
  }
  // Recover the padding to zero
  pad2zero( ptr_necrosis, len );
  
  // 10-3: Find enhancing tumor core
  for( int i = 0; i < len; ++ i ) {
    if( ( ptr_t1ce[ 2 * i ] == T1ce::T1TM && 
        ptr_necrosis[ 2 * i ] != Tumor::NCR ) &&
        ( ptr_t2[ 2 * i ] == T2::T2CSF || 
        ptr_flair[ 2 * i ] == Flair::FTM ) ) {
      ptr_enh[ 2 * i ] = Tumor::ET;
    }
  }
  
  // 10-4: Find edema
  for( int i = 0; i < len; ++ i ) {
    if( ( ptr_flair[ 2 * i ] == Flair::FTM || 
        ptr_t2[ 2 * i ] == T2::T2CSF ||
        ptr_enh[ 2 * i ] == Tumor::ET ) && 
        ( ptr_necrosis[ 2 * i ] != Tumor::NCR &&
        ptr_hemorrhage[ 2 * i ] != Tumor::HMG ) ) {
      ptr_edema[ 2 * i ] = Tumor::ED;
    }
  }
  // 10-5: Find necrosis
  // necrosis enclosed by edema
  inRegion( ptr_enclose_nec, len, ptr_edema, Tumor::ED, 
            ptr_necrosis, Tumor::NCR, 
            region, ptr_nidx, ptr_aidx, nr, nc, ns );
  // Remove necrosis regions separate from edema
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_enclose_nec, ptr_enclose_nec, 
                    1, region ) ) {
      excldRegion( region, ptr_nidx, ptr_enclose_nec, 
                   ptr_edema, Tumor::ED );
    }
  }
  // Recover the padding to zero
  pad2zero( ptr_enclose_nec, len );
  // Extend in ptr_necrosis
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_enclose_nec, ptr_necrosis, 
                    Tumor::NCR, region ) ) {
      extRegion( region, ptr_enclose_nec, 1, .5 );
    }
  }
  // Recover the padding to zero
  pad2zero( ptr_enclose_nec, len );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_enclose_nec[ 2 * i ] == 1 ) {
      ptr_necrosis[ 2 * i ] = Tumor::NCR;
      ptr_edema[ 2 * i ] = 0;
    } else {
      ptr_necrosis[ 2 * i ] = 0;
    }
  }
  
  // 10-6: Find hemorrhage
  // hemorrhage enclosed by edema
  inRegion( ptr_enclose_hem, len, ptr_edema, Tumor::ED, 
            ptr_hemorrhage, Tumor::HMG,
            region, ptr_nidx, ptr_aidx, nr, nc, ns );
  // Remove hemorrhage regions separate from edema
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_enclose_hem, ptr_enclose_hem,
                    1, region ) ) {
      excldRegion( region, ptr_nidx, ptr_enclose_hem,
                   ptr_edema, Tumor::ED );
    }
  }
  // Recover the padding to zero
  pad2zero( ptr_enclose_hem, len );
  // Extend in ptr_hemorrhage
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_enclose_hem, ptr_hemorrhage,
                    Tumor::HMG, region ) ) {
      extRegion( region, ptr_enclose_hem, 1, .5 );
    }
  }
  // Recover the padding to zero
  pad2zero( ptr_enclose_hem, len );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_enclose_hem[ 2 * i ] == 1 ) {
      ptr_hemorrhage[ 2 * i ] = Tumor::HMG;
      ptr_edema[ 2 * i ] = 0;
    } else {
      ptr_hemorrhage[ 2 * i ] = 0;
    }
  }
  pad2zero( ptr_hemorrhage, len );
  // 10-7: HGG or LGG
  int n_enh = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_enh[ 2 * i ] == Tumor::ET ) {
      ++ n_enh;
    }
  }
  if( n_enh > min_enh ) {
    ptr_hgg[ 0 ] = 1;
  } else {
    ptr_hgg[ 0 ] = 0; 
  }
  if( ptr_hgg[ 0 ] == 1 ) {
    // HGG
    // 10-8.1: Find necrosis
    // necrosis enclosed by enh
    inRegion( ptr_enclose_ncr, len, ptr_enh, Tumor::ET, 
              ptr_edema, Tumor::ED, 
              region, ptr_nidx, ptr_aidx, nr, nc, ns );
    for( int i = 0; i < len; ++ i ) {
      if( ptr_enclose_ncr[ 2 * i ] == 1 && 
          ptr_enh[ 2 * i ] != Tumor::ET ) {
        ptr_necrosis[ 2 * i ] = Tumor::NCR;
        ptr_edema[ 2 * i ] = 0;
      }
    }
    // 10-8.2: Remove enh from edema
    for( int i = 0; i < len; ++ i ) {
      if( ptr_enh[ 2 * i ] == Tumor::ET ) {
        ptr_edema[ 2 * i ] = 0;
      }
    }
    // Wrap up the segmentation result
    wrapUp( len, ptr_hemorrhage, ptr_necrosis, ptr_enh, ptr_edema,
            ptr_seg, ptr_tumor );
    // 10-8.3: Remove 3D connected regions with enh.size < min_enh
    for( int i = 0; i < len; ++ i ) {
      if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                      1, region ) ) {
        excldRegion( region, ptr_tumor,
                     ptr_seg, Seg::SET, min_enh );
      }
    }
    pad2zero( ptr_tumor, len );
    for( int i = 0; i < len; ++ i ) {
      if( ptr_tumor[ 2 * i ] == 0 ) {
        ptr_seg[ 2 * i ] = 0;
      }
    }
    
  } else {
    // LGG
    // 10-9.1
    // Wrap up results for LGG
    wrapUp( len, ptr_hemorrhage, ptr_necrosis, ptr_enh,
            ptr_edema, ptr_flair, ptr_seg, ptr_tumor );
    // Remove 3D connected regions with size < min_tumor or 
    // Keep the largest tumor region
    int max_size = 0;
    for( int i = 0; i < len; ++ i ) {
      if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                      1, region ) ) {
        if( region.size() > max_size ) {
          max_size = region.size();
        }
      }
    }
    pad2zero( ptr_tumor, len );
    if( max_size > min_tumor ) {
      for( int i = 0; i < len; ++ i ) {
        if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                        1, region ) ) {
          excldRegion( region, ptr_seg, min_tumor );
        }
      }
    } else {
      int size = max_size - 1;
      for( int i = 0; i < len; ++ i ) {
        if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                        1, region ) ) {
          excldRegion( region, ptr_seg, size );
        }
      }
    }
    furtherSeg( ptr_hgg, len, ptr_seg, 0.50 );
  }
  
  // Restore the segmentation result to a image with the original 
  // dimension
  restoreImg( ptr_idx, ptr_seg, ptr_res_image, len );
  
  SEXP res = PROTECT( allocVector( VECSXP, 2 ) );
  SET_VECTOR_ELT( res, 0, res_image );
  SET_VECTOR_ELT( res, 1, hgg );
  
  SEXP names = PROTECT( allocVector( STRSXP, 2 ) );
  SET_STRING_ELT( names, 0, mkChar( "seg" ) );
  SET_STRING_ELT( names, 1, mkChar( "hgg" ) );
  
  setAttrib( res, R_NamesSymbol, names );
  
  delete [] ptr_hemorrhage;
  delete [] ptr_necrosis;
  delete [] ptr_enh;
  delete [] ptr_edema;
  delete [] ptr_enclose_nec;
  delete [] ptr_enclose_hem;
  delete [] ptr_enclose_ncr;
  delete [] ptr_seg;
  delete [] ptr_tumor;
  
  UNPROTECT( 4 );
  return res;
}