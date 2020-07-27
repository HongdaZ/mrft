#ifndef INITPARMHEALTH_H
#define INITPARMHEALTH_H

#include <vector>
#include <map>
#include <list>

using std::vector;
using std::map;
using std::list;

// Initialize parameters for t1ce and flair images
void initParmHealth3( map<int, vector<double>> &health_parm,
                     int *ptr_seg, const double *ptr_m,
                     const double *ptr_nu2, const double *ptr_intst, 
                     const double *ptr_lambda2, const int *ptr_nidx,
                     const double *ptr_nintst, const double *ptr_alpha,
                     const double *ptr_beta, int len, 
                     int maxit = 20 );

#endif