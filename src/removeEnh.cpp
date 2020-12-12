#include "removeEnh.h"
#include "cnctRegion.h"
#include "pad2zero.h"
#include "excldRegion.h"
#include "tissueType.h"
#include "spread.h"
#include "radius.h"

void removeEnh( int *ptr_seg, const int &label,
                const double &prop_brain, const double &prop_radius, 
                const double &spread_rm, const int &len,
                vector<int> &region, const int *ptr_nidx, 
                const int *ptr_aidx,
                const int &nr, const int &nc, const int &ns ) {
  int n_enh = 0, last_slice = 0, slice = 0, max_n = 0;
  double bound = 0, r_enh = 0, r_brain = 0,
    spread_enh = 0, spread_brain = 0, max_dist = 0;
  // Find number of voxels in each slice
  vector<int> count( ns, 0 );
  for( int i = 0; i < len; ++ i ) {
    ++ count[ ptr_aidx[ 3 * i + 2 ] - 1 ];
  }
  for( int i = 0; i < ns; ++ i ) {
    if( max_n < count[ i ] ) {
      max_n = count[ i ];
    }
  }
  bound = max_n * prop_brain;
  for( int i = 0; i < ns; ++ i ) {
    if( count[ i ] < bound ) {
      last_slice = i + 1;
    } else {
      break;
    }
  }
  for( int i = 0; i < len; ++ i ) {
    slice = ptr_aidx[ 3 * i + Plane::Axial ];
    if( slice < last_slice ) {
      if( cnctRegion( i + 1, ptr_nidx, ptr_aidx, Plane::Axial, 
                      ptr_seg, ptr_seg, 
                      label, region ) ) {
        n_enh = region.size();
        r_enh = radius2D( n_enh );
        spread_enh = spread( region, Plane::Axial, ptr_aidx );
        if( spread_enh > spread_rm ) {
          r_brain = radius2D( count[ slice - 1 ] );
          max_dist = spread_enh * r_enh;
          if( max_dist > prop_radius * r_brain ) {
            excldRegion( region, ptr_seg, 0 );
          }
        }
      }
    }
  }
  pad2zero( ptr_seg, len );
}