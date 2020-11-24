#include "addInside.h"
#include "inRegion.h"
#include "cnctRegion.h"
#include "pTNbr.h"
#include "radius.h"
#include "excldRegion.h"
#include "pad2zero.h"

void addInside( int *ptr_inside, const int &len,
                int *ptr_seg1, const int &label1, 
                int *ptr_seg2, const int &label2,
                vector<int> &region, 
                const int *ptr_nidx, const int *ptr_aidx,
                const int &nr, const int &nc, const int &ns ) {
  int MAX = nr * nc * ns;
  inRegion2D( ptr_inside, len, ptr_seg1, label1, ptr_seg2, label2,
              region, ptr_nidx, ptr_aidx, nr, nc, ns );
  double p_t_nbr;
  double r;
  for( int i = 0; i < len; ++ i ) {
    if( cnctRegion( i + 1, ptr_nidx, ptr_inside, 
                    ptr_inside, 1, region ) ) {
      p_t_nbr = pTNbr( region, ptr_seg1, label1, ptr_nidx );
      r = radius( region.size() );
      if( 1 / p_t_nbr > r / 3 ) {
        excldRegion( region, ptr_inside, MAX );
      }
    }
  }
  pad2zero( ptr_inside, len );
  return;
}