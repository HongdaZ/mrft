#ifndef UPDATESIGMA_H
#define UPDATESIGMA_H

// update sigma for outliers
void updateSigma( int idx,
                  double mu,
                  const double *ptr_intst,
                  double alphal,
                  double betal,
                  double &sigma2 );

#endif