#ifndef UPDATETS_H
#define UPDATETS_H

#include <set>
#include <vector>

using std::set;
using std::vector;

void updateTS( set<int> &region, int current_label,
               double mu,
               double &sigma2,
               double lambda2,
               const int *ptr_seg,
               const int *ptr_nidx,
               const double *ptr_intst,
               const double *ptr_nintst,
               vector<double> &theta,
               double alphal,
               double betal );
                            

#endif