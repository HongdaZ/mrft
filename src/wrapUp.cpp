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
void wrapUp( const int &len, 
             const int *ptr_hemorrhage,
             const int *ptr_necrosis,
             const int *ptr_enh,
             const int *ptr_edema,
             const int *ptr_flair,
             int *ptr_seg, int *ptr_tumor ) {
  for( int i = 0; i < len; ++ i ) {
    if( ptr_hemorrhage[ 2 * i ] == Tumor::HMG ) {
      ptr_seg[ 2 * i ] = Seg::SNET;
      ptr_tumor[ 2 * i ] = 1;
    } else if( ptr_necrosis[ 2 * i ] == Tumor::NCR ) {
      ptr_seg[ 2 * i ] = Seg::SNET;
      ptr_tumor[ 2 * i ] = 1;
    } else if( ptr_enh[ 2 * i ] == Tumor::ET ) {
      // remove enh from LGG
      ptr_seg[ 2 * i ] = 0;
      ptr_tumor[ 2 * i ] = 0;
    } else if( ptr_edema[ 2 * i ] == Tumor::ED ) {
      ptr_tumor[ 2 * i ] = 1;
      if( ptr_flair[ 2 * i ] == Flair::FWM || 
          ptr_flair[ 2 * i ] == Flair::FGM ) {
        ptr_seg[ 2 * i ] = Seg::SNET;
      } else {
        ptr_seg[ 2 * i ] = Seg::SED;
      }
    }
  }
}