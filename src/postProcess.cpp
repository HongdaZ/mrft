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
#include "zeroVector.h"
#include "onRegion.h"

using std::vector;
using std::list;

// Postprocess the results
extern "C" SEXP postProcess( SEXP post_data, SEXP min_enh, 
                             SEXP max_prop_enh_enc,
                             SEXP min_tumor, SEXP spread,
                             SEXP min_prop_tumor_nbr ) ;
SEXP postProcess( SEXP post_data, SEXP min_enh, SEXP max_prop_enh_enc,
                  SEXP min_tumor, SEXP spread,
                  SEXP min_prop_tumor_nbr ) {
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
  const int MAX = nr * nc * ns;
  
  int *ptr_t1ce = INTEGER( t1ce );
  int *ptr_flair = INTEGER( flair );
  int *ptr_t2 = INTEGER( t2 );
  const int *ptr_idx = INTEGER( idx );
  const int *ptr_nidx = INTEGER( nidx );
  const int *ptr_aidx = INTEGER( aidx );
  
  const int len = length( idx );
  SEXP res_image = PROTECT( alloc3DArray( INTSXP, nr, nc, ns ) );
  SEXP code = PROTECT( ScalarInteger( 0 ) );
  
  // Store the results for each tissue type
  int *ptr_seg = new int[ 2 * len ]();
  int *ptr_code = INTEGER( code );
  int *ptr_res_image = INTEGER( res_image );
  int *ptr_hemorrhage = new int[ 2 * len ]();
  int *ptr_necrosis = new int[ 2 * len ]();
  int *ptr_enh = new int[ 2 * len ]();
  int *ptr_edema = new int[ 2 * len ]();
  // FLAIR(4) & T2(4) \ edema \ necrosis \ enh
  // enclosed by edema
  int *ptr_extra_edema = new int[ 2 * len ]();
  // FLAIR(4) & T2(4) \ enh
  int *ptr_whole = new int[ 2 * len ]();
  // Enh inside FLAIR( 4 )
  int *ptr_enclose_enh = new int[ 2 * len ]();
  // Necrosis inside edema
  int *ptr_enclose_nec = new int[ 2 * len ]();
  // Hemorrhage inside edema
  int *ptr_enclose_hem = new int[ 2 * len ]();
  // Necrosis inside enh
  int *ptr_enclose_ncr = new int[ 2 * len ]();
  int *ptr_tumor = new int[ 2 * len ]();
  int *ptr_on = new int[ 2 * len ]();
  // Store the result of findRegion
  vector<int> region;
  region.reserve( len );
  const int m_enh = INTEGER( min_enh )[ 0 ];
  const double m_prop_enh_enc = REAL( max_prop_enh_enc )[ 0 ];
  const int m_tumor = INTEGER( min_tumor )[ 0 ];
  const double spread_factor = REAL( spread )[ 0 ];
  const double m_p_t_nbr = REAL( min_prop_tumor_nbr )[ 0 ];
  
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
  // Recover the padding to zero (cnctRegion)
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
    if( ptr_t2[ 2 * i ] == T2::T2CSF &&
        ptr_t1ce[ 2 * i ] != T1ce::T1TM &&
        ptr_flair[ 2 * i ] != Flair::FTM &&
        ptr_hemorrhage[ 2 * i ] != Tumor::HMG ) {
      ptr_whole[ 2 * i ] = 1;
    }
  }
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_necrosis, ptr_whole,
                    1, region ) ) {
      extRegion( region, ptr_necrosis, Tumor::NCR, 0 );
    }
  }
  // Recover the padding to zero
  pad2zero( ptr_necrosis, len );
  zeroVector( ptr_whole, len );
  // 10-3: Find rough regions of enhancing tumor core
  for( int i = 0; i < len; ++ i ) {
    if( ptr_t1ce[ 2 * i ] == T1ce::T1TM &&
        ( ptr_t2[ 2 * i ] == T2::T2CSF ||
        ptr_flair[ 2 * i ] == Flair::FTM ) &&
        ptr_necrosis[ 2 * i ] != Tumor::NCR &&
        ptr_hemorrhage[ 2 * i ] != Tumor::HMG ) {
      ptr_enh[ 2 * i ] = Tumor::ET;
    }
  }
  // 10-4: Find rough regions of edema
  for( int i = 0; i < len; ++ i ) {
    if( ( ( 
        ptr_t2[ 2 * i ] == T2::T2CSF &&
        ( ptr_flair[ 2 * i ] == Flair::FTM &&
        ptr_t1ce[ 2 * i ] == T1ce::T1GM ) ) ||
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
      extRegion( region, ptr_enclose_nec, 1, .6, true );
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
      extRegion( region, ptr_enclose_hem, 1, .6, true );
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
  // 10-7.1: FLAIR( 4 ) || T2( 4 ) || T1ce( 2 )
  // inside tumor >> edema 
  // Wrap up the segmentation result
  // Now edema only includes edema
  wrapUp( len, ptr_hemorrhage, ptr_necrosis, ptr_enh, ptr_edema,
          ptr_seg, ptr_tumor );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_flair[ 2 * i ] == Flair::FTM ||
        ptr_t2[ 2 * i ] == T2::T2CSF ) {
      ptr_whole[ 2 * i ] = 1;
    } else {
      ptr_whole[ 2 * i ] = 0;
    }
  }
  inRegion( ptr_extra_edema, len, ptr_tumor, 1,
            ptr_whole, 1,
            region, ptr_nidx, ptr_aidx, nr, nc, ns );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_extra_edema[ 2 * i ] == 1 &&
        ptr_seg[ 2 * i ] == 0 ) {
      ptr_tumor[ 2 * i ] = 1;
      ptr_seg[ 2 * i ] = Seg::SED;
      ptr_edema[ 2 * i ] = Tumor::ED;
    }
  }
  zeroVector( ptr_whole, len );
  // 10-7.2: Remove 3D connected regions with
  // size < min_tumor or Keep the tumor regions with 
  // size > max_size * 0.80
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
  if( max_size > m_tumor ) {
    for( int i = 0; i < len; ++ i ) {
      if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                      1, region ) ) {
        excldRegion( region, ptr_seg,  ptr_tumor, ptr_hemorrhage,
                     ptr_necrosis, ptr_enh, ptr_edema, m_tumor );
      }
    }
  } else {
    int size = max_size * 0.80;
    for( int i = 0; i < len; ++ i ) {
      if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                      1, region ) ) {
        excldRegion( region, ptr_seg,  ptr_tumor, ptr_hemorrhage,
                     ptr_necrosis, ptr_enh, ptr_edema, size );
      }
    }
  }
  pad2zero( ptr_tumor, len );
  // // remove regions with prop of enh too large
  // for( int i = 0; i < len; ++ i ) {
  //   if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
  //                   1, region ) ) {
  //     excldRegion( region, ptr_tumor,
  //                  ptr_seg, Seg::SET, -1, m_prop_enh,
  //                  ptr_hemorrhage, ptr_necrosis,
  //                  ptr_enh, ptr_edema );
  //   }
  // }
  // pad2zero( ptr_tumor, len );
  
  // 10-7.3: Extend edema in ( FLAIR( 4 ) && T1ce( 2 ) ) || 
  //                         ( T2( 4 ) && T1ce( 2 ) ) ||
  //                         ( FLAIR( 4 ) && T2( 4 ) && T1ce( 3 ) ) 
  for( int i = 0; i < len; ++ i ) {
    if(
        ( ( ptr_flair[ 2 * i ] == Flair::FTM &&
          ptr_t1ce[ 2 * i ] ==  T1ce::T1GM ) || 
          ( ptr_t2[ 2 * i ] == T2::T2CSF &&
          ptr_t1ce[ 2 * i ] ==  T1ce::T1GM ) || 
        ( ptr_flair[ 2 * i ] == Flair::FTM &&
          ptr_t2[ 2 * i ] == T2::T2CSF &&
          ptr_t1ce[ 2 * i ] ==  T1ce::T1WM ) ) &&
        ptr_tumor[ 2 * i ] == 0 ) {
      ptr_whole[ 2 * i ] = 1;
    } else {
      ptr_whole[ 2 * i ] = 0;
    }
  }
  onRegion( ptr_on, len, 0.5, ptr_tumor, 1, ptr_whole, 1,
            region, spread_factor, m_p_t_nbr,
            ptr_nidx, ptr_aidx, nr, nc, ns );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_on[ 2 * i ] == 1 &&
        ptr_tumor[ 2 * i ] == 0 ) {
      ptr_tumor[ 2 * i ] = 1;
      ptr_seg[ 2 * i ] = Seg::SED;
    }
  }
  // 10-7.4: Add T1ce(4) inside edema
  inRegion( ptr_enclose_enh, len, ptr_tumor, 1,
            ptr_t1ce, T1ce::T1TM,
            region, ptr_nidx, ptr_aidx, nr, nc, ns );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_enclose_enh[ 2 * i ] == 1 &&
        ptr_seg[ 2 * i ] == 0 ) {
      ptr_tumor[ 2 * i ] = 1;
      ptr_seg[ 2 * i ] = Seg::SET;
      ptr_enh[ 2 * i ] = Tumor::ET;
    }
  }
  // 10-7.5: code
  // 0: HGG (no further seg)
  // 1: LGG (further seg)
  // 2: HGG (further seg)
  int n_enh = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_enh[ 2 * i ] == Tumor::ET ) {
      ++ n_enh;
    }
  }
  if( n_enh < m_enh ) {
    // code
    ptr_code[ 0 ] = 1; // LGG (further seg)
    wrapUp( len, ptr_hemorrhage, ptr_necrosis, ptr_enh, ptr_edema,
            ptr_flair, ptr_seg, ptr_tumor );
  } else {
    // 10-7.6: Voxels inside enh is necrosis
    // Find total number of voxels inside convex hull
    // of enh
    // Find whole = enh complement
    for( int i = 0; i < len; ++ i ) {
      if( ptr_t1ce[ 2 * i ] != 0 &&
          ptr_flair[ 2 * i ] != 0 &&
          ptr_t2[ 2 * i ] != 0 &&
          ptr_enh[ 2 * i ] != Tumor::ET ) {
        ptr_whole[ 2 * i ] = 1;
      } else {
        ptr_whole[ 2 * i ] = 0;
      }
    }
    // necrosis enclosed by enh
    inRegion( ptr_enclose_ncr, len, ptr_enh, Tumor::ET,
              ptr_whole, 1,
              region, ptr_nidx, ptr_aidx, nr, nc, ns );
    int n_other = 0, n_enh = 0;
    for( int i = 0; i < len; ++ i ) {
      if( ptr_enclose_ncr[ 2 * i ] == 1 ) {
        ptr_tumor[ 2 * i ] = 1;
        ptr_seg[ 2 * i ] = Seg::SNET;
        ptr_necrosis[ 2 * i ] = Tumor::NCR;
        ptr_edema[ 2 * i ] = 0;
        ptr_hemorrhage[ 2 * i ] = 0;
        ptr_enh[ 2 * i ] = 0;
        ++ n_other;
      } else if( ptr_enh[ 2 * i ] == Tumor::ET ) {
        ++ n_enh;
      }
    }
    int n_tumor = 0;
    for( int i = 0; i < len; ++ i ) {
      if( ptr_tumor[ 2 * i ] == 1 ) {
        ++ n_tumor;
      }
    }
    // Rprintf( "n_other = %d, n_enh = %d\n", n_other, n_enh );
    if( ( n_enh + n_other ) < m_prop_enh_enc *  n_tumor ) {
      ptr_code[ 0 ] = 2; // HGG (further seg)
    } else {
      ptr_code[ 0 ] = 0; // HGG (no further seg)
    }
  }
  // Restore the segmentation result to a image with the original 
  // dimension
  restoreImg( ptr_idx, ptr_seg, ptr_res_image, len );
  
  SEXP res = PROTECT( allocVector( VECSXP, 2 ) );
  SET_VECTOR_ELT( res, 0, res_image );
  SET_VECTOR_ELT( res, 1, code );
  
  SEXP names = PROTECT( allocVector( STRSXP, 2 ) );
  SET_STRING_ELT( names, 0, mkChar( "image" ) );
  SET_STRING_ELT( names, 1, mkChar( "code" ) );
  
  setAttrib( res, R_NamesSymbol, names );
  
  delete [] ptr_hemorrhage;
  delete [] ptr_necrosis;
  delete [] ptr_enh;
  delete [] ptr_edema;
  delete [] ptr_extra_edema;
  delete [] ptr_whole;
  delete [] ptr_enclose_enh;
  delete [] ptr_enclose_nec;
  delete [] ptr_enclose_hem;
  delete [] ptr_enclose_ncr;
  delete [] ptr_seg;
  delete [] ptr_tumor;
  delete [] ptr_on;
  
  UNPROTECT( 4 );
  return res;
}