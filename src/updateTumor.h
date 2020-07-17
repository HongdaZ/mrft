#ifndef UPDATETUMOR_H
#define UPDATETUMOR_H

#include "updateParm.h"

#include <map>
#include <vector>

using std::vector;
using std::map;

// update parameters for tumor regions
void updateTumor( map<int, vector<double>> &tumor_parm, 
                  const map<int, list<int>> &tumor_regions,
                  const int *ptr_seg, const int &len, const double *ptr_m,
                  const double *ptr_nu2, const double *ptr_intst,
                  const double *ptr_lambda2, const int *ptr_nidx,
                  const double *ptr_nintst, const double *ptr_alpha,
                  const double *ptr_beta, const int &maxit );

#endif