#ifndef INITPARM_H
#define INITPARM_H

#include <vector>
#include <map>
#include <list>

using std::vector;
using std::map;
using std::list;

// Initialize parameters
void initParm( const bool first_run, map<int, vector<double>> &health_parm,
               map<int, vector<double>> &tumor_parm,
               int *ptr_seg, const double *ptr_m,
               const double *ptr_nu2, const double *ptr_intst, 
               const double *ptr_lambda2, const int *ptr_nidx,
               const double *ptr_nintst, const double *ptr_alpha,
               const double *ptr_beta,
               map<int, list<int>> &tumor_regions, 
               const double *ptr_a, const double *ptr_b, int len, 
               int maxit = 20 );

#endif