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
extern "C" SEXP postProcess( SEXP post_data, SEXP min_enh, 
                             SEXP max_prop_enh,
                             SEXP min_tumor, SEXP min_prop_net );
SEXP postProcess( SEXP post_data, SEXP min_enh, SEXP max_prop_enh,
                  SEXP min_tumor, SEXP min_prop_net ) {
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
  int *ptr_tmp = new int[ 2 * len ]();
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
  // Store the result of findRegion
  vector<int> region;
  region.reserve( len );
  const int m_enh = INTEGER( min_enh )[ 0 ];
  const double m_prop_enh = REAL( max_prop_enh )[ 0 ];
  const int m_tumor = INTEGER( min_tumor )[ 0 ];
  const double m_prop_net = REAL( min_prop_net )[ 0 ];
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
    if( ptr_t2[ 2 * i ] == T2::T2CSF &&
        ptr_t1ce[ 2 * i ] != T1ce::T1TM &&
        ptr_flair[ 2 * i ] != Flair::FTM &&
        ptr_hemorrhage[ 2 * i ] != Tumor::HMG ) {
      ptr_tmp[ 2 * i ] = 1;
    }
  }
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_necrosis, ptr_tmp,
                    1, region ) ) {
      extRegion( region, ptr_necrosis, Tumor::NCR, 0 );
    }
  }
  // Recover the padding to zero
  pad2zero( ptr_necrosis, len );
  
  // 10-3: Find enhancing tumor core
  // enh enclosed by FLAIR(4)
  // inRegion( ptr_enclose_enh, len, ptr_flair, Flair::FTM,
  //           ptr_t1ce, T1ce::T1TM,
  //           region, ptr_nidx, ptr_aidx, nr, nc, ns );
  // // Extend enh to 3D connected regions
  // for( int i = 0; i < len; ++ i ) {
  //   if( cnctRegion( i + 1, ptr_nidx, ptr_enclose_enh, ptr_t1ce,
  //                   T1ce::T1TM, region ) ) {
  //     excldVoxel( region, ptr_necrosis, Tumor::NCR );
  //     extRegion( region, ptr_enclose_enh, 1, 0.5, true );
  //   }
  // }
  // pad2zero( ptr_enclose_enh, len );
  // for( int i = 0; i < len; ++ i ) {
  //   if( ptr_enclose_enh[ 2 * i ] == 1 ) {
  //     ptr_enh[ 2 * i ] = Tumor::ET;
  //   }
  // }
  for( int i = 0; i < len; ++ i ) {
    if( ptr_t1ce[ 2 * i ] == T1ce::T1TM &&
        ptr_t2[ 2 * i ] == T2::T2CSF &&
        ptr_necrosis[ 2 * i ] != Tumor::NCR &&
        ptr_hemorrhage[ 2 * i ] != Tumor::HMG ) {
      ptr_enh[ 2 * i ] = Tumor::ET;
      
    }
  }
  // 10-4: Find edema
  for( int i = 0; i < len; ++ i ) {
    if( ( ( 
        ptr_flair[ 2 * i ] == Flair::FTM &&
        ptr_t2[ 2 * i ] == T2::T2CSF &&
        ptr_t1ce[ 2 * i ] == T1ce::T1GM ) ||
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
      extRegion( region, ptr_enclose_nec, 1, .5, true );
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
  // 10-6.2: Add T1ce(4) inside Edema
  inRegion( ptr_enclose_enh, len, ptr_edema, Tumor::ED, 
            ptr_t1ce, T1ce::T1TM,
            region, ptr_nidx, ptr_aidx, nr, nc, ns );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_enclose_enh[ 2 * i ] == 1 ) {
      ptr_enh[ 2 * i ] = Tumor::ET;
    }
  }
  // 10-7: HGG or LGG
  int n_enh = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_enh[ 2 * i ] == Tumor::ET ) {
      ++ n_enh;
    }
  }
  if( n_enh > m_enh ) {
    ptr_hgg[ 0 ] = 1;
  } else {
    ptr_hgg[ 0 ] = 0; 
  }
  if( ptr_hgg[ 0 ] == 1 ) {
    // HGG
    // 10-8.1: Find whole = enh complement
    for( int i = 0; i < len; ++ i ) {
      if( ptr_t1ce[ 2 * i ] != 0 &&
          ptr_flair[ 2 * i ] != 0 && 
          ptr_t2[ 2 * i ] != 0 && 
          ptr_enh[ 2 * i ] != Tumor::ET ) {
        ptr_whole[ 2 * i ] = 1;
      }
    }
    // 10-8.2: Find necrosis
    // necrosis enclosed by enh
    inRegion( ptr_enclose_ncr, len, ptr_enh, Tumor::ET, 
              ptr_whole, 1, 
              region, ptr_nidx, ptr_aidx, nr, nc, ns );
    for( int i = 0; i < len; ++ i ) {
      if( ptr_enclose_ncr[ 2 * i ] == 1 && 
          ptr_enh[ 2 * i ] != Tumor::ET ) {
        ptr_necrosis[ 2 * i ] = Tumor::NCR;
        ptr_edema[ 2 * i ] = 0;
      }
    }
    // 10-8.3: Remove enh from edema
    for( int i = 0; i < len; ++ i ) {
      if( ptr_enh[ 2 * i ] == Tumor::ET ) {
        ptr_edema[ 2 * i ] = 0;
      }
    }
    // Wrap up the segmentation result
    wrapUp( len, ptr_hemorrhage, ptr_necrosis, ptr_enh, ptr_edema,
            ptr_seg, ptr_tumor );
    // 10-8.4: Remove 3D connected regions with enh.size < min_enh
    // and percentage of enh > 80%
    for( int i = 0; i < len; ++ i ) {
      if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                      1, region ) ) {
        excldRegion( region, ptr_tumor,
                     ptr_seg, Seg::SET, m_enh, m_prop_enh,
                     ptr_hemorrhage, ptr_necrosis,
                     ptr_enh, ptr_edema );
      }
    }
    pad2zero( ptr_tumor, len );
    // 10-8.5: Find extra edema
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
        ptr_seg[ 2 * i ] = Seg::SED;
        ptr_tumor[ 2 * i ] = 1;
        ptr_edema[ 2 * i ] = Tumor::ED;
      }
    }
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
    if( max_size > m_tumor ) {
      for( int i = 0; i < len; ++ i ) {
        if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                        1, region ) ) {
          excldRegion( region, ptr_seg, m_tumor );
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
    if( max_size > m_tumor ) {
      for( int i = 0; i < len; ++ i ) {
        if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                        1, region ) ) {
          excldRegion( region, ptr_seg, m_tumor );
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
    furtherSeg( ptr_hgg, len, ptr_seg, m_prop_net );
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
  delete [] ptr_tmp;
  delete [] ptr_extra_edema;
  delete [] ptr_whole;
  delete [] ptr_enclose_enh;
  delete [] ptr_enclose_nec;
  delete [] ptr_enclose_hem;
  delete [] ptr_enclose_ncr;
  delete [] ptr_seg;
  delete [] ptr_tumor;
  
  UNPROTECT( 4 );
  return res;
}