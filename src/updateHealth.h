#ifndef UPDATEHEALTH_H
#define UPDATEHEALTH_H

#include "updateParm.h"

#include <map>
#include <vector>

using std::vector;
using std::map;

// debug useage of vector
void updateHealth( map<int, vector<double>> &health_parm,
                   const int *ptr_seg, const int &len, const double *ptr_m,
                   const double *ptr_nu2, const double *ptr_intst,
                   const double *ptr_lambda2, const int *ptr_nidx,
                   const double *ptr_nintst, const double *ptr_alpha,
                   const double *ptr_beta, const int &maxit );

#endif