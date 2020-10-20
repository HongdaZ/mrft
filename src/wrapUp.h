#ifndef WRAPUP_H
#define WRAPUP_H

// Wrap up results for HGG
void wrapUp( const int &len, 
              const int *ptr_hemorrhage,
              const int *ptr_necrosis,
              const int *ptr_enh,
              const int *ptr_edema,
              int *ptr_seg,
              int *ptr_tumor );
// Wrap up results for LGG
void wrapUp( const int &len, int *ptr_hgg, 
             const int *ptr_hemorrhage,
             const int *ptr_necrosis,
             const int *ptr_enh,
             const int *ptr_edema,
             const int *ptr_flair,
             int *ptr_seg );

#endif