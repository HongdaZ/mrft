#include <R.h>
#include <Rinternals.h>

#include "cmpET.h"
#include "nbrLabel.h"
#include "findOutLabel.h"
#include "newTumorLabel.h"

// compare energy for training
void cmpET( int idx, int sc,
            list<map<int, int >> &regions,
            map<int, set<int>> &tumor_regions, 
            set<int> &tumor_labels, set<int> &outl_labels,
            map<int, vector<double>> &health_parm,
            map<int, vector<double>> &tumor_parm, 
            map<int, vector<double>> &outl_parm,
            int *ptr_seg, const int *ptr_nidx,
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
    // New outlier label
    int out_label = findOutLabel( idx, ptr_seg, outl_labels );
    
    // outlier parameters
    double out_mu = - 1, out_sigma2 = 1;
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
    // update parameters for whole region
    if( combine_nrg < split_nrg && sc == 1 ) {
      vector<double> label_whole_parm = region_parm.front();
      int whole_label = label_whole_parm[ 0 ];
      vector<double> whole_parm( ++ label_whole_parm.begin(),
                                 label_whole_parm.end() );
      tumor_parm[ whole_label ] = whole_parm;
    // remove old whole region and add new splitted regions and new outlier
    } else if ( combine_nrg > split_nrg && sc == 1 ) {
      vector<double> label_whole_parm = region_parm.front();
      int whole_label = label_whole_parm[ 0 ];
      
      tumor_labels.erase( whole_label );
      tumor_regions.erase( whole_label );
      tumor_parm.erase( whole_label );
      
      list<vector<double>>::iterator it = ++ region_parm.begin();
      list<map<int, int >>::iterator it_region = ++ regions.begin();
      for( ; it != region_parm.end(); ++ it, ++ it_region ) {
        vector<double>::iterator it_parm = it->begin();
        int new_label = *it_parm;
        vector<double> new_parm( ++ it_parm, it->end() );
        set<int> new_region;
        for( map<int, int>::iterator it_map = it_region->begin();
             it_map != it_region->end(); ++ it_map ) {
          new_region.insert( it_map->first );
        }
        tumor_labels.insert( new_label );
        tumor_regions[ new_label ] = new_region;
        tumor_parm[ new_label ] = new_parm;
      }
      ptr_seg[ 2 * ( idx - 1 ) ] = out_label;
      outl_labels.insert( out_label );
      vector<double> new_outl_parm( ++ outlier_parm.begin(),
                                    outlier_parm.end() );
      outl_parm[ out_label ] = new_outl_parm;
      
      // remove old subregions, parameters and labels,
      // add new whole regions label,region and parameters,
      // set seg to tumor,
      // remove outlier from outl_labels
      // remove outl_parm
    } else if( combine_nrg < split_nrg && sc == 2 ) {
      for( list<vector<double>>::iterator it_list =  ++ region_parm.begin();
           it_list != region_parm.end(); ++ it_list ) {
        int sub_label = it_list->front();
        tumor_labels.erase( sub_label );
        tumor_regions.erase( sub_label );
        tumor_parm.erase( sub_label );
      }
      
      vector<double> label_whole_parm = region_parm.front();
      int whole_label = label_whole_parm[ 0 ];
      set<int> whole_region;
      for( map<int, int>::iterator it_map = regions.begin()->begin();
           it_map != regions.begin()->end(); ++ it_map ) {
        whole_region.insert( it_map->first );
      }
      tumor_labels.insert( whole_label );
      tumor_regions[ whole_label ] = whole_region;
      vector<double> new_whole_parm( ++ region_parm.begin()->begin(),
                                     region_parm.begin()->end() );
      tumor_parm[ whole_label ] = new_whole_parm;
      
      ptr_seg[ 2 * ( idx - 1 ) ] = whole_label;
      outl_labels.erase( out_label );
      outl_parm.erase( out_label );
    // update parameters for subregions
    } else if( combine_nrg > split_nrg && sc == 2 ) {
      for( list<vector<double>>::iterator it_list =  ++ region_parm.begin();
           it_list != region_parm.end(); ++ it_list ) {
        int sub_label = it_list->front();
        vector<double> new_parm( ++ it_list->begin(), it_list->end() );
        tumor_parm[ sub_label ] = new_parm;
      }
    }
  // no split or combine
  } else {
    int curr_label = ptr_seg[ 2 * ( idx - 1 ) ];
    if( curr_label > - 4 && curr_label < 1 ) {
      double energy;
      double min_energy;
      int min_label;
      vector<double> parm;
      double mu;
      double sigma2;
      vector<double> theta;
      for( int i = - 1; i > - 4; -- i ) {
        theta.empty();
        parm = health_parm[ i ]; 
        mu = parm[ 0 ];
        sigma2 = parm[ 1 ];
        theta.insert( theta.begin(), parm.begin() + 2, parm.end() );
        energy = energyY( i, idx, mu, sigma2, ptr_seg, ptr_nidx,
                              ptr_intst, ptr_nintst, theta );
        energy += energyX( i, idx, false, ptr_seg, ptr_nidx, ptr_delta[ 0 ],
                           ptr_gamma[ 0 ] );
        if( i == - 1 ) {
          min_energy = energy;
          min_label = i;
        } else {
          if( energy < min_energy ) {
            min_energy = energy;
            min_label = i;
          }
        } 
      }
      ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
    // tumor or outlier
    } else {
      vector<int> nbr_label = nbrLabel( idx, ptr_seg, ptr_nidx );
      bool have_tumor = false;
      int t_label = 0;
      for( int i = 0; i < 6; ++ i ) {
        if( nbr_label[ i ] != NA_INTEGER && nbr_label[ i ] < - 3 ) {
          have_tumor = true;
          t_label = nbr_label[ i ];
        }
      }
      if( have_tumor ) {
        // tumor energy
        vector<double> t_parm = tumor_parm[ t_label ];
        double t_mu = t_parm[ 0 ];
        double t_sigma2 = t_parm[ 1 ];
        vector<double> t_theta( t_parm.begin() + 2, t_parm.end() );
        double energy_t = energyY( t_label, idx, t_mu, ptr_m[ 2 ], t_sigma2, 
                                   ptr_lambda2[ 3 ], ptr_seg, ptr_nidx,
                                   ptr_intst, ptr_nintst, t_theta, 
                                   ptr_alpha[ 3 ], ptr_beta[ 3 ], ptr_a[ 0 ],
                                   ptr_b[ 0 ] );
        energy_t += energyX( t_label, idx, false, ptr_seg, ptr_nidx,
                             ptr_delta[ 0 ], ptr_gamma[ 0 ] );
        // Outlier energy
        // New outlier label
        int out_label = findOutLabel( idx, ptr_seg, outl_labels );
        
        // outlier parameters
        double out_mu = - 1, out_sigma2 = 1;
        map<int, int> out_region;
        out_region[ idx ] = out_label;
        updateParm( out_mu, out_sigma2, out_region, ptr_m[ 3 ], ptr_m[ 2 ], 
                    ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst, out_label,
                    ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, ptr_alpha[ 3 ],
                    ptr_beta[ 3 ], 20 );
        
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
        vector<double> new_out_parm;
        new_out_parm.push_back( out_mu );
        new_out_parm.push_back( out_sigma2 );
        if( curr_label > 0 && energy_t < out_energy ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = t_label;
          tumor_regions[ t_label ].insert( idx );
          outl_labels.erase( curr_label );
          outl_parm.erase( curr_label );
        } else if( curr_label < - 3 && energy_t > out_energy ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = out_label;
          tumor_regions[ t_label ].erase( idx );
          outl_labels.insert( out_label );
          
          outl_parm[ out_label ] = new_out_parm;
        }
      // not have tumor  
      } else {
        if( curr_label > 0 ) {
          // Find a new tumor region label;
          int new_label = newTumorLabel( 1, tumor_labels );
          // new tumor region parameters
          double mu = -1, sigma2 = 1; // sigma2 has to be non-zero;
          vector<double> theta;
          for( int i = 0; i < 6; ++ i ) {
            theta.push_back( 0 );
          }
          map<int, int> new_region;
          new_region[ idx ] = new_label;
          
          updateParm( mu, theta, sigma2, new_region, ptr_m[ 3 ], ptr_m[ 2 ], 
                      ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst, new_label, 
                      ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, ptr_nintst, 
                      ptr_alpha[ 3 ], ptr_beta[ 3 ], 20 );
          vector<double> new_parm;
          new_parm.push_back( mu );
          new_parm.push_back( sigma2 );
          
          new_parm.insert( new_parm.end(), theta.begin(), theta.end() );
          tumor_parm[ new_label ] = new_parm;

          ptr_seg[ 2 * ( idx - 1 ) ] = new_label;
          tumor_regions[ new_label ].insert( idx ); 
          
          
          outl_labels.erase( curr_label );
          outl_parm.erase( curr_label );
        }
      }
    }
  }
}