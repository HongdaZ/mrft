
#include "updateTumor.h"

// update parameters for tumor regions
void updateTumor( list<int> &tumor_labels,
                  map<int, vector<double>> &tumor_parm, 
                  map<int, list<int>> &tumor_regions,
                  const double *ptr_a, const double *ptr_b,
                  const int *ptr_seg, const double *ptr_m,
                  const double *ptr_intst, const double *ptr_lambda2,
                  const int *ptr_nidx, const double *ptr_nintst,
                  const double *ptr_alpha, const double *ptr_beta, 
                  const int &maxit ) {
  list<int>::iterator it_labels = tumor_labels.begin();
  for( ; it_labels != tumor_labels.end(); ++ it_labels ) {
    list<int> &region = tumor_regions[ *it_labels ];
    if( region.size() > 1 ) {
      // new tumor region parameters
      double t_mu = -1, t_sigma2 = 1; // t_sigma2 has to be non-zero;
      vector<double> t_theta( 6, 0 );
      updateParm( t_mu, t_theta, t_sigma2, region, ptr_m[ 3 ],
                  ptr_m[ 2 ], ptr_a[ 0 ], ptr_b[ 0 ],
                  ptr_intst, *it_labels, ptr_lambda2[ 3 ], ptr_seg,
                  ptr_nidx, ptr_nintst, ptr_alpha[ 3 ],
                  ptr_beta[ 3 ], maxit );
      vector<double> new_parm( 8, 0 );
      new_parm[ 0 ] = t_mu;
      new_parm[ 1 ] = t_sigma2;
      for( int i = 0; i < 6; ++ i ) {
        new_parm[ i + 2 ] = t_theta[ i ];
      }
      tumor_parm[ *it_labels ] = new_parm;
    }
  }
}