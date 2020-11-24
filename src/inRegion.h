#ifndef INREGION_H
#define INREGION_H

#include <vector>

using std::vector;

void inRegion( int *ptr_enclose, const int &len,
               int *ptr_seg1, const int &label1, 
               int *ptr_seg2, const int &label2,
               vector<int> &region, 
               const int *ptr_nidx, const int *ptr_aidx,
               const int &nr, const int &nc, const int &ns );
void inRegion2D( int *ptr_enclose, const int &len,
                 int *ptr_seg1, const int &label1, 
                 int *ptr_seg2, const int &label2,
                 vector<int> &region, 
                 const int *ptr_nidx, const int *ptr_aidx,
                 const int &nr, const int &nc, const int &ns );

#endif