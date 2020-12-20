#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <list>

#include "remove.h"
#include "pad2zero.h"
#include "cnctRegion.h"
#include "excldRegion.h"
#include "tissueType.h"
#include "zeroVector.h"
#include "spread.h"
#include "inRegion.h"
#include "removeSmall.h"

using std::max;
using std::min;
using std::list;

void remove( vector<int> &region, const int *ptr_aidx, 
             const int *ptr_nidx, int *ptr_tumor, int *ptr_seg,
             int *ptr_hemorrhage, int *ptr_necrosis, 
             int *ptr_enh, int *ptr_edema, const int &m_tumor, 
             const int &m_enh, const int &m_enh_enc,
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
  int MAX = nr * nc * ns;
  int n_enc_t = 0;
  int max_size = 0;
  double max_spread_idx = 0, min_spread_idx = 0;
  int n_enh = 0;
  int idx = 0;
  double spread_idx = 0;
  int size = 0;
  double s_factor;
  list<int> region_size;
  list<double> region_spread;

  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                    1, region ) ) {
      spread_idx = spread( region, ptr_aidx );
      region_size.push_back( region.size() );
      if( max_size < region.size() ) {
        max_size = region.size(); 
      }
      region_spread.push_back( spread_idx );
      // min_spread is also updated later
      if( max_spread_idx < spread_idx ) {
        max_spread_idx = spread_idx;
      }
    }
  }
  pad2zero( ptr_tumor, len );
  size = min( m_tumor, max_size );
  min_spread_idx = max_spread_idx;
  
  list<int>::const_iterator it_size = region_size.begin();
  list<double>::const_iterator it_spread = region_spread.begin();
  
  for( ; it_size != region_size.end(); ++ it_size, ++ it_spread ) {
    if( *it_size > 5000 ) {
      spread_idx = *it_spread;
      if( spread_idx < min_spread_idx ) {
        min_spread_idx = spread_idx;
      }
    }
  }
  
  s_factor = max( min_spread_idx, spread_factor );
  // Rprintf( "s_factor = %f\n", s_factor );
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                    1, region ) ) {
      spread_idx = spread( region, ptr_aidx );
      // Rprintf( "spread_idx = %f\n", spread_idx );
      if( region.size() < size &&
          spread_idx > s_factor ) {
        excldRegion( region, ptr_seg,  ptr_tumor, ptr_hemorrhage,
                     ptr_necrosis, ptr_enh, ptr_edema, MAX );
      }
    }
  }
  pad2zero( ptr_tumor, len );
  
  
  // Rprintf( "s_factor = %f\n", s_factor );
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                    1, region ) ) {
      if( region.size() < size ) {
        // Find enh inside the tumor region
        zeroVector( ptr_r_enh, len );
        n_enh = 0;
        for( int j = 0; j < region.size(); ++ j ) {
          idx = region[ j ];
          if( idx != 0 ) {
            if( ptr_enh[ 2 * ( idx - 1 ) ] == Tumor::ET ) {
              ++ n_enh;
              ptr_r_enh[ 2 * ( idx - 1 ) ] = 1;
            }
          }
        }
        if( n_enh < m_enh ) {
          // If not including enh, exclude small sized regions
          if( region.size() < size ) {
            excldRegion( region, ptr_seg,  ptr_tumor, ptr_hemorrhage,
                         ptr_necrosis, ptr_enh, ptr_edema, MAX );
          }
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
                         ptr_necrosis, ptr_enh, ptr_edema, MAX );
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
// Remove small regions
void remove( vector<int> &region, const int *ptr_nidx, int *ptr_tumor,
             int *ptr_seg, int *ptr_hemorrhage, int *ptr_necrosis, 
             int *ptr_enh, int *ptr_edema, const int &m_tumor, 
             const int len ) {
  int max_size = 0;
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_tumor, ptr_tumor,
                    1, region ) ) {
      if( max_size < region.size() ) {
        max_size = region.size(); 
      }
    }
  }
  pad2zero( ptr_tumor, len );
  max_size = min( m_tumor, max_size );
  removeSmall( region, ptr_nidx, ptr_seg, ptr_tumor, ptr_hemorrhage, 
               ptr_necrosis, ptr_enh, ptr_edema,  max_size, len );
}