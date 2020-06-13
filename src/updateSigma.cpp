#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include <set>
#include <vector>

#include "updateSigma.h"

using std::set;
using std::vector;



// region starts from 1
// updateTheta for health regions
void updateSigma( int idx,
                  double mu,
                  const double *ptr_intst,
                  double alphal,
                  double betal,
                  double &sigma2 ) {
  double intst = ptr_intst[ idx - 1 ];
  Rprintf( "sigma2 = %f\n", pow( intst - mu, 2 ) );
  sigma2 = ( betal + pow( intst - mu, 2 ) / 2 ) / ( 1 / 2 + alphal + 1 );
  Rprintf( "sigma2 = %f", sigma2 );
  
  return;
  
}
