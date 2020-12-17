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
#include "remove.h"
#include "assignCSF.h"
#include "removeSmall.h"
#include "addInside.h"
#include "removeSlice.h"
#include "trim.h"
#include "removeEnh.h"
#include "removeEnhBlock.h"
#include "pTNbr.h"

using std::vector;
using std::list;

// Postprocess the results
extern "C" SEXP postProcess( SEXP post_data, SEXP min_enh,
                            SEXP min_enh_enc,
                            SEXP max_prop_enh_enc,
                            SEXP max_prop_enh_slice,
                            SEXP min_tumor, SEXP spread_add,
                            SEXP spread_rm, SEXP spread_trim ) ;
SEXP postProcess( SEXP post_data, SEXP min_enh,
                  SEXP min_enh_enc, SEXP max_prop_enh_enc,
                  SEXP max_prop_enh_slice,
                  SEXP min_tumor, SEXP spread_add,
                  SEXP spread_rm, SEXP spread_trim ) {
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
  SEXP res_edema = PROTECT( alloc3DArray( INTSXP, nr, nc, ns ) );
  SEXP res_csf = PROTECT( alloc3DArray( INTSXP, nr, nc, ns ) );
  list<int> edema_codes, csf_codes;
  int n_edema = 0, n_csf = 0;
  
  // Store the results for each tissue type
  int *ptr_seg = new int[ 2 * len ]();
  int *ptr_edema_regions = new int[ 2 * len ]();
  int *ptr_csf_regions = new int[ 2 * len ]();
  int e_code;
  int *ptr_res_image = INTEGER( res_image );
  int *ptr_res_edema = INTEGER( res_edema );
  int *ptr_res_csf = INTEGER( res_csf );
  int *ptr_hemorrhage = new int[ 2 * len ]();
  int *ptr_necrosis = new int[ 2 * len ]();
  int *ptr_csf = new int[ 2 * len ]();
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
  int *ptr_enclose_csf = new int[ 2 * len ]();
  int *ptr_exclude = new int[ 2 * len ]();
  int *ptr_enh_include = new int[ 2 * len]();
  int *ptr_enh_exclude = new int[ 2 * len]();
  // Each seperate tumor region
  int *ptr_sub_region = new int[ 2 * len ]();
  // copies of ptr_tumor, ptr_seg
  int *ptr_tumor_copy = new int[ 2 * len ]();
  int *ptr_seg_copy = new int[ 2 * len ]();
  int *ptr_local_enh = new int[ 2 * len ]();
  // Store the result of findRegion
  vector<int> region;
  region.reserve( len );
  vector<int> region_tmp;
  region_tmp.reserve( len );
  
  const int m_enh = INTEGER( min_enh )[ 0 ];
  const int m_enh_enc = INTEGER( min_enh_enc )[ 0 ];
  const double m_prop_enh_enc = REAL( max_prop_enh_enc )[ 0 ];
  const double m_prop_enh_slice = REAL( max_prop_enh_slice )[ 0 ];
  const int m_tumor = INTEGER( min_tumor )[ 0 ];
  const double s_add = REAL( spread_add )[ 0 ];
  const double s_rm = REAL( spread_rm )[ 0 ];
  const double s_trim = REAL( spread_trim )[ 0 ];
  
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
    if( ptr_t1ce[ 2 * i ] == T1ce::T1CSF || 
        ptr_flair[ 2 * i ] == Flair::FCSF ) {
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
  for( int i = 0; i < 2 * len; ++ i ) {
    ptr_whole[ i ] = ptr_t1ce[ i ];
  }
  // Remove large 2D slices of T1ce(4)
  removeSlice( ptr_whole, T1ce::T1TM, 20, m_prop_enh_slice, len,
               region, ptr_nidx, ptr_aidx, nr, nc, ns );
  // Remove big regions of enh
  for( int i = 0; i < len; ++ i ) {
    if( ptr_t1ce[ 2 * i ] == T1ce::T1TM ) {
      ptr_enh_include[ 2 * i ] = 1;
    }
  }
  removeEnh( ptr_enh_include, 1, .70, 1 / 2, 3, len, 
             region, ptr_nidx, ptr_aidx, nr, nc, ns );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_whole[ 2 * i ] == T1ce::T1TM &&
        ( ptr_t2[ 2 * i ] == T2::T2CSF ||
        ptr_flair[ 2 * i ] == Flair::FTM ) &&
        ptr_necrosis[ 2 * i ] != Tumor::NCR &&
        ptr_hemorrhage[ 2 * i ] != Tumor::HMG &&
        ptr_enh_include[ 2 * i ] == 1 ) {
      ptr_enh[ 2 * i ] = Tumor::ET;
    }
  }
  zeroVector( ptr_whole, len );
  for( int i = 0; i < len; ++ i ) {
    if( ptr_necrosis[ 2 * i ] == Tumor::NCR ||
        ptr_hemorrhage[ 2 * i ] == Tumor::HMG ) {
      ptr_enh_include[ 2 * i ] = 0;
    }
  }
  // 10-4: Find rough regions of edema
  for( int i = 0; i < 2 * len; ++ i ) {
    ptr_whole[ i ] = ptr_t2[ i ];
    ptr_on[ i ] = ptr_flair[ i ];
  }
  // Remove large 2D slices of T2(4) and Flair(4)
  removeSlice( ptr_whole, T2::T2CSF, 10, 0.8, len,
               region, ptr_nidx, ptr_aidx, nr, nc, ns );
  removeSlice( ptr_on, Flair::FTM, 10, 0.4, len,
               region, ptr_nidx, ptr_aidx, nr, nc, ns );
  for( int i = 0; i < len; ++ i ) {
    if( ( (
        ptr_whole[ 2 * i ] == T2::T2CSF &&
          ( ptr_on[ 2 * i ] == Flair::FTM &&
          ptr_t1ce[ 2 * i ] == T1ce::T1GM ) ) ||
          ptr_enh[ 2 * i ] == Tumor::ET ) &&
          ( ptr_necrosis[ 2 * i ] != Tumor::NCR &&
          ptr_hemorrhage[ 2 * i ] != Tumor::HMG ) ) {
      ptr_edema[ 2 * i ] = Tumor::ED;
    }
  }
  zeroVector( ptr_whole, len );
  zeroVector( ptr_on, len );
  // 10-5: Find necrosis
  // necrosis enclosed by edema
  inRegion2D( ptr_enclose_nec, len, ptr_edema, Tumor::ED,
            ptr_necrosis, Tumor::NCR,
            region, ptr_nidx, ptr_aidx, nr, nc, ns, 1 );
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
    if( cnctRegion( i + 1, ptr_nidx, ptr_aidx, Plane::Axial,
                    ptr_enclose_nec, ptr_necrosis,
                    Tumor::NCR, region ) ) {
      extRegion( region, ptr_enclose_nec, 1, 0.4, false );
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
      extRegion( region, ptr_enclose_hem, 1, .6, false);
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

  // // 10-7.1: Remove 3D connected regions with
  // // size < min_tumor or Keep the tumor regions with
  // // size = max_size
  // wrapUp( len, ptr_hemorrhage, ptr_necrosis, ptr_enh, ptr_edema,
  //         ptr_seg, ptr_tumor );
  // remove( region, ptr_aidx, ptr_nidx, ptr_tumor, ptr_seg,
  //         ptr_hemorrhage, ptr_necrosis, ptr_enh, ptr_edema,
  //         m_tumor, m_enh, m_enh_enc, len, nr, nc, ns, s_rm );
  // // Trim tumor region
  // trim( ptr_tumor, ptr_exclude,
  //       ptr_nidx, ptr_aidx, region, len, s_trim );
  // for( int i = 0; i < len; ++ i ) {
  //   if( ptr_tumor[ 2 * i ] == 0 ) {
  //     ptr_seg[ 2 * i ] = 0;
  //     ptr_hemorrhage[ 2 * i ] = 0;
  //     ptr_necrosis[ 2 * i ] = 0;
  //     ptr_enh[ 2 * i ] = 0;
  //     ptr_edema[ 2 * i ] = 0;
  //   }
  // }
  // removeSmall( region, ptr_nidx, ptr_seg, ptr_tumor,
  //              ptr_hemorrhage, ptr_necrosis, ptr_enh, ptr_edema,
  //              200, len );
  // // 10-7.2: FLAIR( 4 ) || T2( 4 )
  // // inside tumor >> edema
  // // Wrap up the segmentation result
  // // Now edema only includes edema
  // for( int i = 0; i < len; ++ i ) {
  //   if( (
  //       // ( ptr_flair[ 2 * i ] == Flair::FTM &&
  //       // ptr_t1ce[ 2 * i ] ==  T1ce::T1GM ) ||
  //       ( ptr_t2[ 2 * i ] == T2::T2CSF &&
  //         ptr_t1ce[ 2 * i ] ==  T1ce::T1GM ) ||
  //         ( ptr_flair[ 2 * i ] == Flair::FTM &&
  //         ptr_t2[ 2 * i ] == T2::T2CSF &&
  //         ptr_t1ce[ 2 * i ] ==  T1ce::T1WM ) ) ) {
  //     ptr_whole[ 2 * i ] = 1;
  //   } else {
  //     ptr_whole[ 2 * i ] = 0;
  //   }
  // }
  // inRegion( ptr_extra_edema, len, ptr_tumor, 1,
  //           ptr_whole, 1,
  //           region, ptr_nidx, ptr_aidx, nr, nc, ns );
  // for( int i = 0; i < len; ++ i ) {
  //   if( ptr_extra_edema[ 2 * i ] == 1 &&
  //       ptr_seg[ 2 * i ] == 0 ) {
  //     ptr_tumor[ 2 * i ] = 1;
  //     ptr_seg[ 2 * i ] = Seg::SED;
  //     ptr_edema[ 2 * i ] = Tumor::ED;
  //   }
  // }
  // zeroVector( ptr_whole, len );
  // // Remove 3D connected regions with
  // // size < 200
  // removeSmall( region, ptr_nidx, ptr_seg, ptr_tumor,
  //              ptr_hemorrhage, ptr_necrosis, ptr_enh, ptr_edema,
  //              200, len );
  // // 10-7.3: Extend edema in ( FLAIR( 4 ) && T1ce( 2 ) )
  // for( int i = 0; i < len; ++ i ) {
  //   if(
  //     ( ptr_flair[ 2 * i ] == Flair::FTM &&
  //       ptr_t1ce[ 2 * i ] ==  T1ce::T1GM ) &&
  //       ptr_exclude[ 2 * i ] != 1 &&
  //       ptr_tumor[ 2 * i ] == 0 ) {
  //     ptr_whole[ 2 * i ] = 1;
  //   } else {
  //     ptr_whole[ 2 * i ] = 0;
  //   }
  // }
  // onRegion( ptr_on, len, 1, ptr_tumor, 1, ptr_whole, 1,
  //           region, s_add,
  //           ptr_nidx, ptr_aidx, nr, nc, ns );
  // for( int i = 0; i < len; ++ i ) {
  //   if( ptr_on[ 2 * i ] == 1 &&
  //       ptr_tumor[ 2 * i ] == 0 ) {
  //     ptr_tumor[ 2 * i ] = 1;
  //     ptr_seg[ 2 * i ] = Seg::SED;
  //     ptr_edema[ 2 * i ] = Tumor::ED;
  //   }
  // }
  // // ( T2( 4 ) && T1ce( 2 ) )
  // for( int i = 0; i < len; ++ i ) {
  //   if(
  //     ( ptr_t2[ 2 * i ] == T2::T2CSF &&
  //       ptr_t1ce[ 2 * i ] ==  T1ce::T1GM ) &&
  //       ptr_exclude[ 2 * i ] != 1 &&
  //       ptr_tumor[ 2 * i ] == 0 ) {
  //     ptr_whole[ 2 * i ] = 1;
  //   } else {
  //     ptr_whole[ 2 * i ] = 0;
  //   }
  // }
  // onRegion( ptr_on, len, 0.2, ptr_tumor, 1, ptr_whole, 1,
  //           region, s_add,
  //           ptr_nidx, ptr_aidx, nr, nc, ns );
  // for( int i = 0; i < len; ++ i ) {
  //   if( ptr_on[ 2 * i ] == 1 &&
  //       ptr_tumor[ 2 * i ] == 0 ) {
  //     ptr_tumor[ 2 * i ] = 1;
  //     ptr_seg[ 2 * i ] = Seg::SED;
  //     ptr_edema[ 2 * i ] = Tumor::ED;
  //   }
  // }
  // // ( FLAIR( 4 ) && T2( 4 ) && T1ce( 3 ) )
  // for( int i = 0; i < len; ++ i ) {
  //   if(
  //     ( ptr_flair[ 2 * i ] == Flair::FTM &&
  //       ptr_t2[ 2 * i ] == T2::T2CSF &&
  //       ptr_t1ce[ 2 * i ] ==  T1ce::T1WM ) &&
  //       ptr_exclude[ 2 * i ] != 1 &&
  //       ptr_tumor[ 2 * i ] == 0 ) {
  //     ptr_whole[ 2 * i ] = 1;
  //   } else {
  //     ptr_whole[ 2 * i ] = 0;
  //   }
  // }
  // onRegion( ptr_on, len, 0.2, ptr_tumor, 1, ptr_whole, 1,
  //           region, s_add,
  //           ptr_nidx, ptr_aidx, nr, nc, ns );
  // for( int i = 0; i < len; ++ i ) {
  //   if( ptr_on[ 2 * i ] == 1 &&
  //       ptr_tumor[ 2 * i ] == 0 ) {
  //     ptr_tumor[ 2 * i ] = 1;
  //     ptr_seg[ 2 * i ] = Seg::SED;
  //     ptr_edema[ 2 * i ] = Tumor::ED;
  //   }
  // }
  // // Remove large blocks of enh
  // zeroVector( ptr_whole, len ); // tumor \ enh
  // zeroVector( ptr_on, len );    // enh
  // for( int i = 0; i < len; ++ i ) {
  //   if( ptr_tumor[ 2 * i ] == 1 &&
  //       ptr_enh[ 2 * i ] == 0 ) {
  //     ptr_whole[ 2 * i ] = 1;
  //   }
  // }
  // for( int i = 0; i < len; ++ i ) {
  //   if( ptr_enh[ 2 * i ] == Tumor::ET ) {
  //     ptr_on[ 2 * i ] = 1;
  //   }
  // }
  // removeEnhBlock( ptr_enh_exclude, ptr_on, 1,
  //                 ptr_whole, 1, 0.3, 0.8, len, region,
  //                 ptr_nidx, ptr_aidx, nr, nc, ns );
  // for( int i = 0; i < len; ++ i ) {
  //   if( ptr_enh_exclude[ 2 * i ] == 1 ) {
  //     ptr_tumor[ 2 * i ] = 0;
  //     ptr_enh[ 2 * i ] = 0;
  //     ptr_seg[ 2 * i ] = 0;
  //   }
  // }
  // zeroVector( ptr_whole, len );
  // zeroVector( ptr_on, len );
  // // T1ce(4) || ( FLAIR( 4 ) && T2( 4 ) && T1ce( 3 ) )
  // for( int i = 0; i < len; ++ i ) {
  //   if( ( ptr_enh_include[ 2 * i ] == 1 ||
  //         ( ptr_flair[ 2 * i ] == Flair::FTM &&
  //           ptr_t2[ 2 * i ] == T2::T2CSF &&
  //           ptr_t1ce[ 2 * i ] ==  T1ce::T1WM ) ) &&
  //       ptr_exclude[ 2 * i ] != 1 &&
  //       ptr_tumor[ 2 * i ] == 0 ) {
  //     ptr_whole[ 2 * i ] = 1;
  //   } else {
  //     ptr_whole[ 2 * i ] = 0;
  //   }
  // }
  // onRegion( ptr_on, len, 0.5, ptr_tumor, 1, ptr_whole, 1,
  //           region, 5,
  //           ptr_nidx, ptr_aidx, nr, nc, ns );
  // for( int i = 0; i < len; ++ i ) {
  //   if( ptr_on[ 2 * i ] == 1 &&
  //       ptr_tumor[ 2 * i ] == 0 ) {
  //     ptr_tumor[ 2 * i ] = 1;
  //     if( ptr_enh_include[ 2 * i ] == 1 ) {
  //       ptr_seg[ 2 * i ] = Seg::SET;
  //       ptr_enh[ 2 * i ] = Tumor::ET;
  //     } else {
  //       ptr_seg[ 2 * i ] = Seg::SED;
  //       ptr_edema[ 2 * i ] = Tumor::ED;
  //     }
  //   }
  // }
  // zeroVector( ptr_whole, len );
  // zeroVector( ptr_on, len );
  // // 10-7.3.2 Add voxels inside tumor (2D)
  // zeroVector( ptr_whole, len );
  // zeroVector( ptr_extra_edema, len );
  // 
  // for( int i = 0; i < len; ++ i ) {
  //   if( ptr_tumor[ 2 * i ] == 0 ) {
  //     ptr_whole[ 2 * i ] = 1;
  //   }
  // }
  // addInside( ptr_extra_edema, len, ptr_tumor, 1, ptr_whole, 1,
  //            region, ptr_nidx, ptr_aidx, nr, nc, ns );
  // for( int i = 0; i < len; ++ i ) {
  //   if( ptr_extra_edema[ 2 * i ] == 1 &&
  //       ptr_tumor[ 2 * i ] == 0 ) {
  //     ptr_tumor[ 2 * i ] = 1;
  //     if( ptr_t1ce[ 2 * i ] == T1ce::T1TM
  //         //   &&
  //         // ( ptr_flair[ 2 * i ] == Flair::FTM ||
  //         // ptr_t2[ 2 * i ] == T2::T2CSF )
  //         ) {
  //       ptr_seg[ 2 * i ] = Seg::SET;
  //       ptr_enh[ 2 * i ] = Tumor::ET;
  //     } else {
  //       ptr_seg[ 2 * i ] = Seg::SED;
  //       ptr_edema[ 2 * i ] = Tumor::ED;
  //     }
  //   }
  // }
  // zeroVector( ptr_whole, len );
  // zeroVector( ptr_extra_edema, len );
  // // 10-7.3.3: T1ce(3) in tumor and connected to enh is enh
  // for( int i = 0; i < len; ++ i ) {
  //   if( ptr_tumor[ 2 * i ] == 1 &&
  //       ptr_t1ce[ 2 * i ] == T1ce::T1WM ) {
  //     ptr_whole[ 2 * i ] = 1;
  //   }
  // }
  // for( int i = 0; i < len; ++ i ) {
  //   if( cnctRegion( i + 1, ptr_nidx, ptr_whole, ptr_whole,
  //                   1, region ) ) {
  //     excldRegion( region, ptr_nidx, ptr_whole,
  //                  ptr_seg, Seg::SET );
  //   }
  // }
  // for( int i = 0; i < len; ++ i ) {
  //   if( ptr_whole[ 2 * i ] == 1 ) {
  //     ptr_enh[ 2 * i ] = Tumor::ET;
  //     ptr_hemorrhage[ 2 * i ] = 0;
  //     ptr_necrosis[ 2 * i ] = 0;
  //     ptr_edema[ 2 * i ] = 0;
  //     ptr_seg[ 2 * i ] = Seg::SET;
  //   }
  // }
  // zeroVector( ptr_whole, len );
  // // 10-7.4: code
  // // 0: HGG (no further seg)
  // // 1: LGG (further seg)
  // // 2: HGG (further seg)
  // for( int i = 0; i < len; ++ i ) {
  //   ptr_tumor_copy[ 2 * i ] = ptr_tumor[ 2 * i ];
  // }
  // for( int j = 0; j < len; ++ j ) {
  //   if( cnctRegion( j + 1, ptr_nidx, ptr_tumor_copy, ptr_tumor_copy, 1,
  //                   region_tmp ) ) {
  //     int idx = 0;
  //     zeroVector( ptr_sub_region, len );
  //     for( vector<int>::const_iterator it = region_tmp.begin();
  //          it != region_tmp.end(); ++ it ) {
  //       idx = *it;
  //       if( idx != 0 ) {
  //         ptr_sub_region[ 2 * ( idx - 1 ) ] = 1;
  //       }
  //     }
  //     vector<int> local_enh ( region_tmp.begin(), region_tmp.end() );
  //     int n_enh = 0;
  //     for( vector<int>::iterator it = local_enh.begin();
  //          it != local_enh.end(); ++ it ) {
  //       idx = *it;
  //       if( idx != 0 ) {
  //         if( ptr_enh[ 2 * ( idx - 1 ) ] == 0 ) {
  //           *it = 0;
  //         } else {
  //           ++ n_enh;
  //         }
  //       }
  //     }
  //     double p_edema_nbr = pTNbr( local_enh, ptr_edema,
  //                                 Tumor::ED, ptr_nidx );
  //     double p_nt_nbr = pTNbr( local_enh, ptr_tumor, 0, ptr_nidx );
  //     if( n_enh < m_enh &&
  //         p_nt_nbr < p_edema_nbr ) {
  //       // code
  //       e_code = 1; // LGG (further seg)
  //     } else {
  //       // 10-7.5: tumor || T1ce( 1 ) || Flair( 4 )
  //       // || T2( 4 ) \ T1ce( 4 )
  //       // inside enh is necrosis
  //       // Find total number of voxels inside convex hull
  //       // of enh
  //       // Find whole = enh complement
  //       zeroVector( ptr_local_enh, len );
  //       for( vector<int>::iterator it = local_enh.begin();
  //            it != local_enh.end(); ++ it ) {
  //         idx = *it;
  //         if( idx != 0 ) {
  //           ptr_local_enh[ 2 * ( idx - 1 ) ] = Tumor::ET;
  //         }
  //       }
  //       zeroVector( ptr_whole, len );
  //       for( int i = 0; i < len; ++ i ) {
  //         if( ( ptr_tumor[ 2 * i ] == 1 ||
  //             ptr_t1ce[ 2 * i ] == T1ce::T1CSF ||
  //             ptr_flair[ 2 * i ] == Flair::FTM ||
  //             ptr_t2[ 2 * i ] == T2::T2CSF ) &&
  //             ptr_enh[ 2 * i ] != Tumor::ET ) {
  //           ptr_whole[ 2 * i ] = 1;
  //         } else {
  //           ptr_whole[ 2 * i ] = 0;
  //         }
  //       }
  // 
  //       // necrosis enclosed by enh
  //       zeroVector( ptr_enclose_ncr, len );
  //       inRegion( ptr_enclose_ncr, len, ptr_local_enh, Tumor::ET,
  //                 ptr_whole, 1,
  //                 region, ptr_nidx, ptr_aidx, nr, nc, ns );
  //       // Remove new regions with size > 200
  //       zeroVector( ptr_on, len );
  //       for( int i = 0; i < len; ++ i ) {
  //         if( ptr_tumor[ 2 * i ] == 0 &&
  //             ptr_enclose_ncr[ 2 * i ] == 1 ) {
  //           ptr_on[ 2 * i ] = 1;
  //         } else {
  //           ptr_on[ 2 * i ] = 0;
  //         }
  //       }
  //       for( int i = 0; i < len; ++ i ) {
  //         if( cnctRegion( i + 1, ptr_nidx, ptr_enclose_ncr, ptr_on,
  //                         1, region ) ) {
  //           excldRegion( region, ptr_enclose_ncr, 200 );
  //         }
  //       }
  //       pad2zero( ptr_enclose_csf, len );
  //       int n_other = 0, n_enh = 0;
  //       for( int i = 0; i < len; ++ i ) {
  //         if( ptr_enclose_ncr[ 2 * i ] == 1 ) {
  //           ptr_sub_region[ 2 * i ] = 1;
  //           ptr_tumor[ 2 * i ] = 1;
  //           ptr_seg[ 2 * i ] = Seg::SNET;
  //           ptr_necrosis[ 2 * i ] = Tumor::NCR;
  //           ptr_edema[ 2 * i ] = 0;
  //           ptr_hemorrhage[ 2 * i ] = 0;
  //           ptr_enh[ 2 * i ] = 0;
  //           ++ n_other;
  //         } else if( ptr_local_enh[ 2 * i ] == Tumor::ET ) {
  //           ++ n_enh;
  //         }
  //       }
  //       int n_tumor = 0;
  //       for( int i = 0; i < len; ++ i ) {
  //         if( ptr_sub_region[ 2 * i ] == 1 &&
  //             ptr_tumor[ 2 * i ] == 1 ) {
  //           ++ n_tumor;
  //         }
  //       }
  //       // Rprintf( "n_other = %d, n_enh = %d\n", n_other, n_enh );
  //       if( ( n_enh + n_other ) < m_prop_enh_enc *  n_tumor ||
  //           p_nt_nbr > p_edema_nbr ) {
  //         e_code = 2; // HGG (further seg)
  //       } else {
  //         e_code = 3; // HGG (no further seg)
  //       }
  //     }
  //     for( int i = 0; i < len; ++ i ) {
  //       if( ptr_sub_region[ 2 * i ] == 1 ) {
  //         ptr_edema_regions[ 2 * i ] = e_code;
  //       }
  //     }
  //   }
  // }
  // pad2zero( ptr_tumor_copy, len );
  // 
  // // 10-7.6: Add T1ce(4) inside edema
  // assignCSF( ptr_csf, ptr_t1ce, ptr_flair, ptr_tumor, len );
  // zeroVector( ptr_whole, len );
  // for( int i = 0; i < len; ++ i ) {
  //   if( ptr_t1ce[ 2 * i ] == T1ce::T1TM &&
  //       ptr_csf[ 2 * i ] == 0 ) {
  //     ptr_whole[ 2 * i ] = 1;
  //   } else {
  //     ptr_whole[ 2 * i ] = 0;
  //   }
  // }
  // for( int i = 0; i < len; ++ i ) {
  //   ptr_tumor_copy[ 2 * i ] = ptr_tumor[ 2 * i ];
  // }
  // for( int j = 0; j < len; ++ j ) {
  //   if( cnctRegion( j + 1, ptr_nidx, ptr_tumor_copy, ptr_tumor_copy, 1,
  //                   region_tmp ) ) {
  //     e_code = ptr_edema_regions[ 2 * j ];
  //     int idx = 0;
  //     zeroVector( ptr_sub_region, len );
  //     for( vector<int>::const_iterator it = region_tmp.begin();
  //          it != region_tmp.end(); ++ it ) {
  //       idx = *it;
  //       if( idx != 0 ) {
  //         ptr_sub_region[ 2 * ( idx - 1 ) ] = 1;
  //       }
  //     }
  //     zeroVector( ptr_enclose_enh, len );
  //     inRegion( ptr_enclose_enh, len, ptr_sub_region, 1,
  //               ptr_whole, 1,
  //               region, ptr_nidx, ptr_aidx, nr, nc, ns );
  //     // Remove new enh regions separate from old enh
  //     for( int i = 0; i < len; ++ i ) {
  //       if( cnctRegion( i + 1, ptr_nidx, ptr_aidx, Plane::Axial,
  //                       ptr_enclose_enh,
  //                       ptr_enclose_enh, 1, region ) ) {
  //         excldRegion( region, ptr_nidx, ptr_enclose_enh,
  //                      ptr_enh, Tumor::ET );
  //       }
  //     }
  //     pad2zero( ptr_enclose_enh, len );
  //     // // Extend in ptr_enh_include
  //     // for( int i = 0; i < len; ++ i ) {
  //     //   if( cnctRegion( i + 1, ptr_nidx, ptr_enclose_enh,
  //     //                   ptr_enh_include, 1, region ) ) {
  //     //     extRegion( region, ptr_enclose_enh, 1, .6, false);
  //     //   }
  //     // }
  //     // // Recover the padding to zero
  //     // pad2zero( ptr_enclose_enh, len );
  //     for( int i = 0; i < len; ++ i ) {
  //       if( ptr_enclose_enh[ 2 * i ] == 1 &&
  //           ptr_seg[ 2 * i ] == 0 ) {
  //         ptr_sub_region[ 2 * i ] = 1;
  //         ptr_tumor[ 2 * i ] = 1;
  //         ptr_seg[ 2 * i ] = Seg::SET;
  //         ptr_enh[ 2 * i ] = Tumor::ET;
  //       }
  //     }
  //     for( int i = 0; i < len; ++ i ) {
  //       if( ptr_sub_region[ 2 * i ] == 1 ) {
  //         ptr_edema_regions[ 2 * i ] = e_code;
  //       }
  //     }
  //   }
  // }
  // pad2zero( ptr_tumor_copy, len );
  // removeSmall( region, ptr_nidx, ptr_seg, ptr_tumor,
  //              ptr_hemorrhage, ptr_necrosis, ptr_enh, ptr_edema,
  //              200, len );
  wrapUp( len, ptr_edema_regions,
          ptr_hemorrhage, ptr_necrosis, ptr_enh, ptr_edema,
          ptr_flair, ptr_seg, ptr_tumor );

  
  // Restore the segmentation result to a image with the original
  // dimension
  SEXP edema_code = PROTECT( allocVector( INTSXP, n_edema )  );
  int *ptr_edema_code_ = INTEGER( edema_code );
  list<int>::const_iterator it_e = edema_codes.begin();
  
  for( int i = 0; it_e != edema_codes.end(); ++ it_e, ++ i ) {
    ptr_edema_code_[ i ] = *it_e;
  }
  SEXP csf_code = PROTECT( allocVector( INTSXP, n_csf )  );
  int *ptr_csf_code_ = INTEGER( csf_code );
  list<int>::const_iterator it_c = csf_codes.begin();
  
  for( int i = 0; it_c != csf_codes.end(); ++ it_c, ++ i ) {
    ptr_csf_code_[ i ] = *it_c;
  }
  
  restoreImg( ptr_idx, ptr_seg, ptr_res_image, len );
  restoreImg( ptr_idx, ptr_edema_regions, ptr_res_edema, len );
  restoreImg( ptr_idx, ptr_csf_regions, ptr_res_csf, len );
  SEXP res = PROTECT( allocVector( VECSXP, 5 ) );
  SET_VECTOR_ELT( res, 0, res_image );
  SET_VECTOR_ELT( res, 1, res_edema );
  SET_VECTOR_ELT( res, 2, edema_code );
  SET_VECTOR_ELT( res, 3, res_csf );
  SET_VECTOR_ELT( res, 4, csf_code );
  
  SEXP names = PROTECT( allocVector( STRSXP, 5 ) );
  SET_STRING_ELT( names, 0, mkChar( "image" ) );
  SET_STRING_ELT( names, 1, mkChar( "edema" ) );
  SET_STRING_ELT( names, 2, mkChar( "edema_code" ) );
  SET_STRING_ELT( names, 3, mkChar( "csf" ) );
  SET_STRING_ELT( names, 4, mkChar( "csf_code" ) );
  
  setAttrib( res, R_NamesSymbol, names );
  
  delete [] ptr_hemorrhage;
  delete [] ptr_necrosis;
  delete [] ptr_csf;
  delete [] ptr_enh;
  delete [] ptr_edema;
  delete [] ptr_extra_edema;
  delete [] ptr_whole;
  delete [] ptr_enclose_enh;
  delete [] ptr_enclose_nec;
  delete [] ptr_enclose_hem;
  delete [] ptr_enclose_ncr;
  delete [] ptr_seg;
  delete [] ptr_edema_regions;
  delete [] ptr_csf_regions;
  delete [] ptr_tumor;
  delete [] ptr_on;
  delete [] ptr_enclose_csf;
  delete [] ptr_exclude;
  delete [] ptr_enh_include;
  delete [] ptr_enh_exclude;
  delete [] ptr_sub_region;
  delete [] ptr_tumor_copy;
  delete [] ptr_seg_copy;
  delete [] ptr_local_enh;
  
  UNPROTECT( 7 );
  return res;
}