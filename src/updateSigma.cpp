#include <R.h>
#include <Rinternals.h>

#include "updateSigma.h"

void updateSigma( int idx,
                  double mu,
                  const double *ptr_intst,
                  double alphal,
                  double betal,
                  double &sigma2 ) {
  double intst = ptr_intst[ idx - 1 ];
  sigma2 = ( betal + pow( intst - mu, 2 ) / 2 ) / ( .5 + alphal + 1 );
  return;
  
}
