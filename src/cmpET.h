#ifndef CMPET_H
#define CMPET_H


#include <vector>
#include <map>
#include <list>

#include "updateParm.h"
#include "energyY.h"
#include "energyX.h"

using std::vector;
using std::map;
using std::list;

// compare energy for training
void cmpET( vector<int> &region, int idx, int sc,
            const vector<int> &labels, const vector<int> &regions,
            vector<int> &tumor_labels, vector<int> &outl_labels,
            const vector<double> &health_parm,
            vector<double> &tumor_parm, 
            vector<double> &outl_parm,
            int *ptr_seg, const int *ptr_nidx,
            const double *ptr_intst, const double *ptr_nintst,
            const double *ptr_delta, const double *ptr_gamma, 
            const double *ptr_alpha, const double *ptr_beta, 
            const double *ptr_lambda2, const double *ptr_a,
            const double *ptr_b, const double *ptr_m,
            const double *ptr_nu2,
            // outlier_parm( 3, 0 )
            vector<double> &outlier_parm,
            // theta( 6, 0 )
            vector<double> &theta,
            // tmp_parm( 9, 0 )
            vector<double> &tmp_parm,
            // out_theta( 6, 0 )
            vector<double> &out_theta,
            // new_out_parm( 2, 0 )
            vector<double> &new_out_parm,
            // whole_parm( 8, 0 )
            vector<double> &whole_parm,
            int &n_tumor, int &n_outl,
            const int &len );

#endif