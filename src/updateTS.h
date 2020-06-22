#ifndef UPDATETS_H
#define UPDATETS_H

#include <set>
#include <vector>

using std::set;
using std::vector;

// region starts from 1
// updateTheta and sigma2 for health and tumorous regions
void updateTS( const set<int> &region, int curr_label,
               const double mu,
               double &sigma2,
               const double lambda2,
               const int *ptr_seg,
               const int *ptr_nidx,
               const double *ptr_intst,
               const double *ptr_nintst,
               vector<double> &theta,
               const double alphal,
               const double betal );
                            

#endif