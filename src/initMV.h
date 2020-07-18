#ifndef INITMV_H
#define INITMV_H

#include <list>

using std::list;

// initialize yln_, yln_i and yl_
void initMV( const list<int> &region, double *yln_, double *yln_i,
                  double *yl_, double &sum_y, const double *ptr_intst, 
                  const int *ptr_nidx, const double *ptr_nintst, 
                  const int *ptr_seg, const int curr_label );
// Initialize yln and yl;
void initMV( double *yl, const double *yl_, double *yln, const double *yln_,
             const double *yln_i, const int nrow, const double mu );


#endif