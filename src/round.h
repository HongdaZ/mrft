#ifndef ROUND_H
#define ROUND_H

// get roundness of the intersection between seg1 and seg2
double round( const int *ptr_seg1, const int &label1, 
              const int *ptr_seg2, const int &label2,
              const int &len, const int *ptr_nidx, 
              const int *ptr_aidx );
#endif