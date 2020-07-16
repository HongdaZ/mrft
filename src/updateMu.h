#ifndef UPDATEMU_H
#define UPDATEMU_H

#include <list>
#include <vector>
#include <map>

using std::list;
using std::vector;
using std::map;

// update mu for healthy cells
double updateMu( const list<int> &region,  
                 const double sigma2,
                 const double m,
                 const double nu2,
                 const vector<double> theta,
                 const double *ptr_intst );
// update mu for tumor cells
double updateMu( const list<int> &region,  
                 const double sigma2,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const vector<double> theta,
                 const double *ptr_intst );
// update mu for outliers
double updateMu( const list<int> &region,  
                 const double sigma2,
                 const double m,
                 const double mk_1,
                 const double a,
                 const double b,
                 const double *ptr_intst );

#endif