#ifndef UPDATEMU_H
#define UPDATEMU_H

#include <set>
#include <vector>
#include <map>

using std::set;
using std::vector;
using std::map;

// update mu for healthy cells
double updateMu( set<int> &region,  
                 double sigma2,
                 double m,
                 double nu2,
                 vector<double> theta,
                 const double *ptr_intst );
// update mu for tumor cells
double updateMu( map<int, int> &region,  
                 double sigma2,
                 double m,
                 double mk,
                 double a,
                 double b,
                 vector<double> theta,
                 const double *ptr_intst );
// update mu for outliers
double updateMu( map<int, int> &region,  
                 double sigma2,
                 double m,
                 double mk_1,
                 double a,
                 double b,
                 const double *ptr_intst );

#endif