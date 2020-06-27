#ifndef CMPET_H
#define CMPET_H


#include <set>
#include <vector>
#include <map>
#include <list>

#include "updateParm.h"
#include "energyY.h"
#include "energyX.h"

using std::set;
using std::vector;
using std::map;
using std::list;

// compare energy for training
void cmpET( int idx, int sc,
            list<map<int, int >> &regions,
            map<int, set<int>> &tumor_regions, 
            set<int> &tumor_labels, set<int> &outl_labels,
            map<int, vector<double>> &health_parm,
            map<int, vector<double>> &tumor_parm, 
            map<int, vector<double>> &outl_parm,
            int *ptr_seg, const int *ptr_nidx,
            const double *ptr_intst, const double *ptr_nintst,
            const double *ptr_delta, const double *ptr_gamma, 
            const double *ptr_alpha, const double *ptr_beta, 
            const double *ptr_lambda2, const double *ptr_a,
            const double *ptr_b, const double *ptr_m,
            const double *ptr_nu2 );

#endif