#ifndef UPDATESIGMA_H
#define UPDATESIGMA_H

#include <set>
#include <vector>

using std::set;
using std::vector;

// update sigma for outliers
void updateSigma( int idx,
                  double mu,
                  const double *ptr_intst,
                  double alphal,
                  double betal,
                  double &sigma2 );

#endif