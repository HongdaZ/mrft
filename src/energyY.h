#ifndef ENERGYY_H
#define ENERGYY_H

#include <list>
#include <vector>

using std::list;
using std::vector;

double energyY( const list<int> &region,
                double mu,
                double mk1,
                double sigma2,
                double lambda2,
                const int *ptr_seg,
                const int *ptr_nidx,
                const double *ptr_intst,
                const double *ptr_nintst,
                const vector<double> &theta,
                double alphak,
                double betak,
                double a, double b );
double energyY( const int current_label,
                const int current_idx,
                double mu,
                double sigma2,
                const int *ptr_seg,
                const int *ptr_nidx,
                const double *ptr_intst,
                const double *ptr_nintst,
                const vector<double> &theta );
// Single point energyY for tumor or outlier
double energyY( const int curr_label,
                const int curr_idx,
                double mu,
                double mk1,
                double sigma2,
                double lambda2,
                const int *ptr_seg,
                const int *ptr_nidx,
                const double *ptr_intst,
                const double *ptr_nintst,
                const vector<double> &theta,
                double alphak,
                double betak,
                double a, double b );
                            

#endif