#ifndef UPDATETUMOR_H
#define UPDATETUMOR_H

#include "updateParm.h"

#include <map>
#include <vector>

using std::vector;
using std::map;

// update parameters for tumor regions
void updateTumor( list<int> &tumor_labels,
                  map<int, vector<double>> &tumor_parm, 
                  map<int, list<int>> &tumor_regions,
                  const double *ptr_a, const double *ptr_b,
                  int *ptr_seg, const double *ptr_m,
                  const double *ptr_intst, const double *ptr_lambda2,
                  const int *ptr_nidx, const double *ptr_nintst,
                  const double *ptr_alpha, const double *ptr_beta, 
                  const int &maxit );

#endif