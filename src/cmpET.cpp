#include <R.h>
#include <Rinternals.h>

#include "cmpET.h"
#include "findOutLabel.h"
#include "newTumorLabel.h"
#include "eraseRegion.h"
#include "eraseOutl.h"
#include "addRegion.h"
#include "addOutl.h"
#include "tumorNbr.h"

// compare energy for training
void cmpET( int idx, int sc,
            list<int> &labels, list<list<int>> &regions,
            map<int, list<int>> &tumor_regions, 
            vector<int> &tumor_labels, vector<int> &outl_labels,
            map<int, vector<double>> &health_parm,
            map<int, vector<double>> &tumor_parm, 
            map<int, vector<double>> &outl_parm,
            int *ptr_seg, const int *ptr_nidx,
            const double *ptr_intst, const double *ptr_nintst,
            const double *ptr_delta, const double *ptr_gamma, 
            const double *ptr_alpha, const double *ptr_beta, 
            const double *ptr_lambda2, const double *ptr_a,
            const double *ptr_b, const double *ptr_m,
            const double *ptr_nu2,
            // outlier_parm( 3, 0 )
            vector<double> &outlier_parm,
            // theta( 6, 0 )
            vector<double> &theta,
            // tmp_parm( 9, 0 )
            vector<double> &tmp_parm,
            // out_theta( 6, 0 )
            vector<double> &out_theta,
            // new_out_parm( 2, 0 )
            vector<double> &new_out_parm,
            // whole_parm( 8, 0 )
            vector<double> &whole_parm,
            int &n_tumor, int &n_outl ) {
  
  if( sc != 0 ) {  // split or combine
    // energy of whole region, subregions and outlier
    int len = regions.size();
    vector<double> nrg( len + 1, 0 ); 
    // parameters for regions
    list<vector<double>> region_parm;
    list<list<int>>::iterator it = regions.begin();
    list<int>::iterator it_label = labels.begin();
    vector<double>::iterator it_nrg = nrg.begin();
    for( ; it != regions.end(); ++ it, ++ it_label, ++ it_nrg ) {
      double mu = - 1, sigma2 = 1;
      int curr_label = *it_label;
      updateParm( mu, theta, sigma2, *it, ptr_m[ 3 ], ptr_m[ 2 ], 
                  ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst, curr_label, 
                  ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, ptr_nintst, 
                  ptr_alpha[ 3 ], ptr_beta[ 3 ], 20 );
      // curr_label, mu, sigma2, theta
      tmp_parm[ 0 ] = curr_label;
      tmp_parm[ 1 ] = mu;
      tmp_parm[ 2] = sigma2;
      for( int i = 0; i < 6; ++ i ) {
        tmp_parm[ i + 3 ] = theta[ i ];
      }
      region_parm.push_back( tmp_parm );
      
      // calculate energy
      double energy = energyY( *it, mu, ptr_m[ 2 ], sigma2, 
                               ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, 
                               ptr_intst, ptr_nintst, 
                               theta, ptr_alpha[ 3 ], ptr_beta[ 3 ], 
                               ptr_a[ 0 ], ptr_b[ 0 ] );
      *it_nrg = energy;
    }
    // New outlier label
    int out_label = findOutLabel( idx, ptr_seg, outl_labels );
    
    // outlier parameters
    double out_mu = - 1, out_sigma2 = 1;
    list<int> out_region;
    out_region.push_back( idx );
    updateParm( out_mu, out_sigma2, out_region, ptr_m[ 3 ], ptr_m[ 2 ], 
                ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst, out_label,
                ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, ptr_alpha[ 3 ],
                ptr_beta[ 3 ], 20 );
    outlier_parm[ 0 ] = out_label;
    outlier_parm[ 1 ] = out_mu;
    outlier_parm[ 2 ] = out_sigma2;
    
    double out_energy = energyX( out_label, idx, true, ptr_seg, ptr_nidx,
                                 ptr_delta[ 0 ], ptr_gamma[ 0 ] );
    vector<double> out_theta( 6, 0 );
    
    out_energy += energyY( out_label, idx, out_mu, ptr_m[ 2 ], out_sigma2,
                           ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, ptr_intst,
                           ptr_nintst, out_theta, ptr_alpha[ 3 ], 
                           ptr_beta[ 3 ], ptr_a[ 0 ], ptr_b[ 0 ] );
    *it_nrg = out_energy;
    it_nrg = nrg.begin();
    double combine_nrg = *it_nrg;
    double split_nrg = 0;
    for( ++ it_nrg; it_nrg != nrg.end(); ++ it_nrg ) {
      split_nrg += *it_nrg;
    }
    vector<double> &label_whole_parm = region_parm.front();
    int whole_label = label_whole_parm[ 0 ];
    // update parameters for whole region
    if( combine_nrg <= split_nrg && sc == 1 ) {
      for( int i = 0; i < 8; ++ i ) {
        whole_parm[ i ] = label_whole_parm[ 1 + i ];
      }
      tumor_parm[ whole_label ] = whole_parm;
    // remove old whole region and add new splitted regions and new outlier
    } else if ( combine_nrg > split_nrg && sc == 1 ) {
      eraseRegion( whole_label, tumor_labels, tumor_regions, tumor_parm, 
                   n_tumor );
      list<vector<double>> new_region_parm( ++ region_parm.begin() ,
                                            region_parm.end() );
      list<list<int>> new_regions( ++ regions.begin(), 
                                        regions.end() );
      addRegion( ptr_seg, new_region_parm, new_regions, tumor_labels, tumor_regions, 
                 tumor_parm, n_tumor );
      
      ptr_seg[ 2 * ( idx - 1 ) ] = out_label;
      
      for( int i = 0; i < 2; ++ i ) {
        new_out_parm[ i ] = outlier_parm[ 1 + i ];
      }
      addOutl( out_label, new_out_parm, outl_labels, outl_parm, n_outl );
      // if( out_label == 32 ) {
      //   Rprintf( "place 1: idx = %d; mu = %f, sigma2 = %f\n", idx, new_out_parm[ 0 ], 
      //            new_out_parm[ 1 ] );
      // }
      // remove old subregions, parameters and labels,
      // add new whole regions label,region and parameters,
      // set seg to tumor,
      // remove outlier from outl_labels
      // remove outl_parm
    } else if( combine_nrg < split_nrg && sc == 2 ) {
      for( list<vector<double>>::iterator it_list =  ++ region_parm.begin();
           it_list != region_parm.end(); ++ it_list ) {
        int sub_label = it_list->front();
        eraseRegion( sub_label, tumor_labels, tumor_regions, tumor_parm, 
                     n_tumor );
      }
      
      list<vector<double>> new_whole_parm;
      list<list<int>> new_whole_region;
      new_whole_parm.push_back( label_whole_parm );
      new_whole_region.push_back( regions.front() );
      addRegion( ptr_seg, new_whole_parm, new_whole_region, tumor_labels, 
                 tumor_regions, tumor_parm, n_tumor );
      ptr_seg[ 2 * ( idx - 1 ) ] = whole_label;
      eraseOutl( out_label, outl_labels, outl_parm, n_outl );
      // if( out_label == 32 ) {
      //   bool in = outl_parm.find( 32 ) != outl_parm.end();
      //   Rprintf( "place 5: idx = %d, in list = %d\n", 
      //            idx, in );
      // }
      
    // update parameters for subregions
    } else if( combine_nrg >= split_nrg && sc == 2 ) {
      for( list<vector<double>>::iterator it_list =  ++ region_parm.begin();
           it_list != region_parm.end(); ++ it_list ) {
        int sub_label = it_list->front();
        vector<double> &new_parm = whole_parm; 
        for( int i = 0; i < 8; ++ i ) {
          new_parm[ i ] = *( it_list->begin() + i + 1); 
        }
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
      double mu;
      double sigma2;
      for( int i = - 1; i > - 4; -- i ) {
        vector<double> &parm = health_parm[ i ]; 
        mu = parm[ 0 ];
        sigma2 = parm[ 1 ];
        for( int i = 0; i < 6; ++ i ) {
          theta[ i ] = parm[ i + 2 ];
        }
        energy = energyY( i, idx, mu, sigma2, ptr_seg, ptr_nidx,
                          ptr_intst, ptr_nintst, theta );
        energy += energyX( i, idx, false, ptr_seg, ptr_nidx, 
                           ptr_delta[ 0 ], ptr_gamma[ 0 ] );
        // Rprintf( "h_label = %d; h_energy = %f; h_mu = %f; h_sigma2 = %f;\n",
        //          i, energy, mu, sigma2 );
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
      // check if having tumor neighbor
      int tumor_idx = firstTumorNbr( idx, ptr_seg, ptr_nidx ) ;
      int t_label;
      vector<double> &t_theta = theta;
      // have tumor neighbor
      if( tumor_idx != 0 ) {
        int t_idx = tumor_idx;
        t_label = ptr_seg[ 2 * ( t_idx - 1 ) ];
        // tumor energy
        vector<double> &t_parm = tumor_parm[ t_label ];
        double t_mu = t_parm[ 0 ];
        double t_sigma2 = t_parm[ 1 ];
        for( int i = 0; i < 6; ++ i ) {
          t_theta[ i ] = t_parm[ i + 2 ];
        }
        double t_energy = energyY( t_label, idx, t_mu, ptr_m[ 2 ], 
                                   t_sigma2, ptr_lambda2[ 3 ], ptr_seg, 
                                   ptr_nidx, ptr_intst, ptr_nintst, 
                                   t_theta, ptr_alpha[ 3 ], ptr_beta[ 3 ], 
                                   ptr_a[ 0 ], ptr_b[ 0 ] );
        t_energy += energyX( t_label, idx, false, ptr_seg, ptr_nidx,
                             ptr_delta[ 0 ], ptr_gamma[ 0 ] );
        // Rprintf( "t_label = %d; t_energy = %f; t_mu = %f; t_sigma2 = %f;\n",
        //          t_label, t_energy, t_mu, t_sigma2 );
        // Outlier energy
        // New outlier label
        int out_label = findOutLabel( idx, ptr_seg, outl_labels );
        
        
        double out_energy = energyX( out_label, idx, false, ptr_seg, 
                                     ptr_nidx, ptr_delta[ 0 ], 
                                     ptr_gamma[ 0 ] );
        // outlier parameters
        double out_mu = - 1, out_sigma2 = 1;
        list<int> out_region;
        out_region.push_back( idx );
        if( curr_label > 0 ) {
          new_out_parm = outl_parm[ curr_label ];
          out_mu = new_out_parm[ 0 ];
          out_sigma2 = new_out_parm[ 1 ];
        } else {
          updateParm( out_mu, out_sigma2, out_region, ptr_m[ 3 ], 
                      ptr_m[ 2 ], ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst, 
                      out_label, ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, 
                      ptr_alpha[ 3 ], ptr_beta[ 3 ], 20 );
          new_out_parm[ 0 ] = out_mu;
          new_out_parm[ 1 ] = out_sigma2;
        }
        
        out_energy += energyY( out_label, idx, out_mu, ptr_m[ 2 ], 
                               out_sigma2, ptr_lambda2[ 3 ], ptr_seg, 
                               ptr_nidx, ptr_intst, ptr_nintst, out_theta, 
                               ptr_alpha[ 3 ],ptr_beta[ 3 ], ptr_a[ 0 ],
                               ptr_b[ 0 ] );
        // Rprintf( "out_label = %d; out_energy = %f; out_mu = %f; out_sigma2 = %f\n",
        //          out_label, out_energy, out_mu, out_sigma2 );
        if( curr_label > 0 && t_energy <= out_energy ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = t_label;
          tumor_regions[ t_label ].remove( idx );
          eraseOutl( curr_label, outl_labels, outl_parm, n_outl );
          // if( curr_label == 32 ) {
          //   bool in = outl_parm.find( 32 ) != outl_parm.end();
          //   Rprintf( "place 3: idx = %d, in list = %d\n", idx, in );
          // }
        } else if( curr_label < - 3 && t_energy >= out_energy ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = out_label;
          tumor_regions[ t_label ].remove( idx );
          addOutl( out_label, new_out_parm, outl_labels, outl_parm, 
                   n_outl );
          // if( out_label == 32 ) {
          //   Rprintf( "place 2: idx = %d; mu = %f, sigma2 = %f\n", idx, new_out_parm[ 0 ], 
          //            new_out_parm[ 1 ] );
          // }
        }
        // not have tumor  
      } else {
        if( curr_label > 0 ) {
          // Find a new tumor region label;
          int new_label = newTumorLabel( 1, tumor_labels );
          // new tumor region parameters
          double mu = -1, sigma2 = 1; // sigma2 has to be non-zero;
          list<int> new_region;
          list<list<int>> new_regions;
          new_region.push_back( idx );
          new_regions.push_back( new_region );
          
          updateParm( mu, theta, sigma2, new_region, ptr_m[ 3 ], 
                      ptr_m[ 2 ], ptr_a[ 0 ], ptr_b[ 0 ], 
                      ptr_intst, new_label, ptr_lambda2[ 3 ], ptr_seg, 
                      ptr_nidx, ptr_nintst, ptr_alpha[ 3 ], 
                      ptr_beta[ 3 ], 20 );
          vector<double> &new_parm = tmp_parm;
          list<vector<double>> new_region_parm;
          new_parm[ 0 ] = new_label;
          new_parm[ 1 ] = mu;
          new_parm[ 2 ] = sigma2;
          for( int i = 0; i < 6; ++ i ) {
            new_parm[ i + 3 ] = theta[ i ];
          }
          new_region_parm.push_back( new_parm );
          
          addRegion( ptr_seg, new_region_parm, new_regions, tumor_labels, 
                     tumor_regions, tumor_parm, n_tumor );
          ptr_seg[ 2 * ( idx - 1 ) ] = new_label;
          eraseOutl( curr_label, outl_labels, outl_parm, n_outl );
          // if( curr_label == 32 ) {
          //   bool in = outl_parm.find( 32 ) != outl_parm.end();
          //   Rprintf( "place 4: idx = %d, in list = %d\n", idx, in );
          // }
        }
      }
    }
  }
}
