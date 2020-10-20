#include "tissueType.h"
#include "wrapUp.h"

// Wrap up results for HGG
void wrapUp( const int &len, 
             const int *ptr_hemorrhage,
             const int *ptr_necrosis,
             const int *ptr_enh,
             const int *ptr_edema,
             int *ptr_seg,
             int *ptr_tumor ) {
  for( int i = 0; i < len; ++ i ) {
    if( ptr_hemorrhage[ 2 * i ] == Tumor::HMG ) {
      ptr_seg[ 2 * i ] = Seg::SNET;
      ptr_tumor[ 2 * i ] = 1;
    } else if( ptr_necrosis[ 2 * i ] == Tumor::NCR ) {
      ptr_seg[ 2 * i ] = Seg::SNET;
      ptr_tumor[ 2 * i ] = 1;
    } else if( ptr_enh[ 2 * i ] == Tumor::ET ) {
      ptr_seg[ 2 * i ] = Seg::SET;
      ptr_tumor[ 2 * i ] = 1;
    } else if( ptr_edema[ 2 * i ] == Tumor::ED ) {
      ptr_seg[ 2 * i ] = Seg::SED;
      ptr_tumor[ 2 * i ] = 1;
    }
  }
}

// Wrap up results for LGG
void wrapUp( const int &len, int *ptr_hgg, 
             const int *ptr_hemorrhage,
             const int *ptr_necrosis,
             const int *ptr_enh,
             const int *ptr_edema,
             const int *ptr_flair,
             int *ptr_seg ) {
  int n_edema = 0, n_flair23 = 0;
  for( int i = 0; i < len; ++ i ) {
    if( ptr_hemorrhage[ 2 * i ] == Tumor::HMG ) {
      ptr_seg[ 2 * i ] = Seg::SNET;
    } else if( ptr_necrosis[ 2 * i ] == Tumor::NCR ) {
      ptr_seg[ 2 * i ] = Seg::SNET;
    } else if( ptr_enh[ 2 * i ] == Tumor::ET ) {
      ptr_seg[ 2 * i ] = Seg::SET;
    } else if( ptr_edema[ 2 * i ] == Tumor::ED ) {
      ++ n_edema;
      if( ptr_flair[ 2 * i ] == Flair::FWM || 
          ptr_flair[ 2 * i ] == Flair::FGM ) {
        ptr_seg[ 2 * i ] = Seg::SNET;
        ++ n_flair23;
      } else {
        ptr_seg[ 2 * i ] = Seg::SED;
      }
    }
  }
  if( n_flair23 > ( 0.25 * n_edema ) ) {
    ptr_hgg[ 0 ] = -1;
  } else {
    ptr_hgg[ 0 ] = -2;
  }
}