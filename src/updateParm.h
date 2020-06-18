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
                 set<int> &region,
                 double m,
                 double nu2,
                 const double *ptr_intst,
                 int curr_label,
                 double lambda2,
                 const int *ptr_seg,
                 const int *ptr_nidx,
                 const double *ptr_nintst,
                 double alphal,
                 double betal,
                 int maxit = 20  );
// update parameters for tumor regions
void updateParm( double &mu, vector<double> &theta, double &sigma2,
                 map<int, int> &region,
                 double m, double mk, double a, double b,
                 const double *ptr_intst,
                 int curr_label,
                 double lambda2,
                 const int *ptr_seg,
                 const int *ptr_nidx,
                 const double *ptr_nintst,
                 double alphal,
                 double betal,
                 int maxit = 20  );
// update parameters for outliers
void updateParm( double &mu, double &sigma2,
                 map<int, int> &region,
                 double m, double mk, double a, double b,
                 const double *ptr_intst,
                 int curr_label,
                 double lambda2,
                 const int *ptr_seg,
                 const int *ptr_nidx,
                 const double *ptr_nintst,
                 double alphal,
                 double betal,
                 int idx,
                 int maxit = 20 );

#endif