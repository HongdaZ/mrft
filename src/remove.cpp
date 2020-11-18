#include "remove.h"
#include "pad2zero.h"
#include "cnctRegion.h"
#include "excldRegion.h"
#include "tissueType.h"
#include "zeroVector.h"
#include "spread.h"
#include "inRegion.h"

void remove( vector<int> &region, const int *ptr_aidx, 
             const int *ptr_nidx, int *ptr_tumor, int *ptr_seg,
             int *ptr_hemorrhage, int *ptr_necrosis, 
             int *ptr_enh, int *ptr_edema, const int &m_tumor, 
             const int &m_enh_enc,
             const int len, const int &nr, const int &nc, const int &ns,
             const double &spread_factor ) {
  int *ptr_r_enh = new int[ 2 * len ]();
  int *ptr_enclose_tumor = new int[ 2 * len ]();
  int *ptr_tmp_tumor = new int[ 2 * len ]();
  vector<int> tmp_region;
  tmp_region.reserve( len );
  
  for( int i = 0; i < len; ++ i ) {
    ptr_tmp_tumor[ 2 * i ] = ptr_tumor[ 2 * i ];
  }
  int n_enc_t = 0;
  int max_size = 0;
  bool have_enh = false;
  int idx = 0;
  double spread_idx = 0;
  int size = 0;
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
    size = m_tumor;
  } else {
    size = max_size;
  }
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                    1, region ) ) {
      // Find enh inside the tumor region
      zeroVector( ptr_r_enh, len );
      have_enh = false;
      for( int j = 0; j < region.size(); ++ j ) {
        idx = region[ j ];
        if( idx != 0 ) {
          if( ptr_enh[ 2 * ( idx - 1 ) ] == Tumor::ET ) {
            have_enh = true;
            ptr_r_enh[ 2 * ( idx - 1 ) ] = 1;
          }
        }
      }
      if( ! have_enh ) {
        // If not including enh, exclude small sized regions
        excldRegion( region, ptr_seg,  ptr_tumor, ptr_hemorrhage,
                     ptr_necrosis, ptr_enh, ptr_edema, size );
      } else {
        // Remove tumor regions with spread > 2
        spread_idx = spread( region, ptr_aidx );
        if( spread_idx > spread_factor ) {
          excldRegion( region, ptr_seg,  ptr_tumor, ptr_hemorrhage,
                       ptr_necrosis, ptr_enh, ptr_edema, size );
        } else {
          // Remove regions with tumor enclosed enh < m_enh_enc
          n_enc_t = 0;
          zeroVector( ptr_enclose_tumor, len );
          inRegion( ptr_enclose_tumor, len, ptr_r_enh, 1,
                    ptr_tmp_tumor, 1, tmp_region, ptr_nidx, ptr_aidx,
                    nr, nc, ns );
          for( int k = 0; k < len; ++ k ) {
            if( ptr_enclose_tumor[ 2 * k ] == 1 ) {
              ++ n_enc_t;
            }
          }
          if( n_enc_t < m_enh_enc ) {
            excldRegion( region, ptr_seg,  ptr_tumor, ptr_hemorrhage,
                         ptr_necrosis, ptr_enh, ptr_edema, size );
          }
        }
      }
    }
  }
  pad2zero( ptr_tumor, len );
  
  delete [] ptr_r_enh;
  delete [] ptr_enclose_tumor;
  delete [] ptr_tmp_tumor;
}