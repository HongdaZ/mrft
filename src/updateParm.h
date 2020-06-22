#ifndef UPDATEPARM_H
#define UPDATEPARM_H


#include <set>
#include <vector>

#include "updateMu.h"
#include "updateTS.h"
#include "updateSigma.h"

using std::set;
using std::vector;

// update parameters for healthy cells
void updateParm( double &mu, vector<double> &theta, double &sigma2, 
                 const set<int> &region,
                 const double m,
                 const double nu2,
                 const double *ptr_intst,
                 int curr_label,
                 const double lambda2,
                 const int *ptr_seg,
                 const int *ptr_nidx,
                 const double *ptr_nintst,
                 const double alphal,
                 const double betal,
                 int maxit );
// update parameters for tumor regions
void updateParm( double &mu, vector<double> &theta, double &sigma2, 
                 const map<int, int> &region,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double *ptr_intst,
                 const int curr_label,
                 const double lambda2,
                 const int *ptr_seg,
                 const int *ptr_nidx,
                 const double *ptr_nintst,
                 const double alphal,
                 const double betal,
                 int maxit );
// update parameters for outliers
void updateParm( double &mu, double &sigma2, 
                 const map<int, int> &region,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double *ptr_intst,
                 int curr_label,
                 const double lambda2,
                 const int *ptr_seg,
                 const int *ptr_nidx,
                 const double alphal,
                 const double betal,
                 int maxit );

#endif