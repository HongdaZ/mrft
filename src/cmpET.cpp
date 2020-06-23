#include <R.h>
#include <Rinternals.h>

#include "cmpET.h"

// compare energy for training
void cmpET( int idx, int sc,
            list<map<int, int >> &regions,
            map<int, set<int>> &tumor_regions, 
            set<int> &tumor_labels, set<int> &outl_labels,
            map<int, vector<double>> &health_parm,
            map<int, vector<double>> &tumor_parm, 
            map<int, vector<double>> &outl_parm,
            int *ptr_seg, const int *ptr_idx, const int *ptr_nidx,
            const double *ptr_intst, const double *ptr_nintst,
            const double *ptr_delta, const double *ptr_gamma, 
            const double *ptr_alpha, const double *ptr_beta, 
            const double *ptr_lambda2, const double *ptr_a,
            const double *ptr_b, const double *ptr_m,
            const double *ptr_nu2 ) {
  
  if( sc != 0 ) {  // split or combine
    // energy of whole region, subregions and outlier
    vector<double> nrg; 
    // parameters for regions
    list<vector<double>> region_parm;
    // parameters for outler
    vector<double> outlier_parm;
    for( list<map<int, int>>::iterator it = regions.begin();
         it != regions.end(); ++ it ) {
      double mu = - 1, sigma2 = 1;
      vector<double> theta;
      for( int i = 0; i < 6; ++ i ) {
        theta.push_back( 0 );
      }
      map<int, int>::iterator it_map = it->begin();
      int curr_label = it_map->second;
      updateParm( mu, theta, sigma2, *it, ptr_m[ 3 ], ptr_m[ 2 ], ptr_a[ 0 ],
                  ptr_b[ 0 ], ptr_intst, curr_label, ptr_lambda2[ 3 ], 
                  ptr_seg, ptr_nidx, ptr_nintst, ptr_alpha[ 3 ], 
                  ptr_beta[ 3 ], 20 );
      // curr_label, mu, sigma2, theta
      vector<double> tmp_parm; 
      tmp_parm.push_back( curr_label );
      tmp_parm.push_back( mu );
      tmp_parm.push_back( sigma2 );
      
      tmp_parm.insert( tmp_parm.end(), theta.begin(), theta.end() );
      region_parm.push_back( tmp_parm );
      
      // calculate energy
      set<int> t_region;
      for( ; it_map != it->end(); ++ it_map ) {
        t_region.insert( it_map->first );
      } 
      double energy = energyY( t_region, mu, ptr_m[ 2 ], sigma2, 
                               ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, ptr_intst,
                               ptr_nintst, theta, ptr_alpha[ 3 ], 
                               ptr_beta[ 3 ], ptr_a[ 0 ], ptr_b[ 0 ] );
      nrg.push_back( energy );
    }
    // outlier parameters
    double out_mu = - 1, out_sigma2 = 1;
    int out_label = *( -- outl_labels.end() ) + 1;
    map<int, int> out_region;
    out_region[ idx ] = out_label;
    updateParm( out_mu, out_sigma2, out_region, ptr_m[ 3 ], ptr_m[ 2 ], 
                ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst, out_label,
                ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, ptr_alpha[ 3 ],
                ptr_beta[ 3 ], 20 );
    outlier_parm.push_back( out_label );
    outlier_parm.push_back( out_mu );
    outlier_parm.push_back( out_sigma2 );
    
    double out_energy = energyX( out_label, idx, true, ptr_seg, ptr_nidx,
                                 ptr_delta[ 0 ], ptr_gamma[ 0 ] );
    vector<double> out_theta;
    for( int i = 0; i < 6; ++ i ) {
      out_theta.push_back( 0 );
    }
    
    out_energy += energyY( out_label, idx, out_mu, ptr_m[ 2 ], out_sigma2,
                           ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, ptr_intst,
                           ptr_nintst, out_theta, ptr_alpha[ 3 ], 
                           ptr_beta[ 3 ], ptr_a[ 0 ], ptr_b[ 0 ] );
    nrg.push_back( out_energy );
    vector<double>::iterator it_nrg = nrg.begin();
    double combine_nrg = *it_nrg;
    double split_nrg = 0;
    for( ++ it_nrg; it_nrg != nrg.end(); ++ it_nrg ) {
      split_nrg += *it_nrg;
    }
    //if()
  }
}