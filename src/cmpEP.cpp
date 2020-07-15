#include <R.h>
#include <Rinternals.h>

#include "cmpEP.h"
#include "nbrLabel.h"
#include "findOutLabel.h"
#include "newTumorLabel.h"
#include "eraseRegion.h"
#include "eraseOutl.h"
#include "addRegion.h"
#include "addOutl.h"

// compare energy for prediction
void cmpEP( int idx, int sc,
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
            const double *ptr_nu2,
            // outlier_parm( 3, 0 )
            vector<double> &outlier_parm,
            // theta( 6, 0 )
            vector<double> &theta,
            // tmp_parm( 9, 0 )
            vector<double> tmp_parm,
            // out_theta( 6, 0 )
            vector<double> &out_theta,
            // new_out_parm( 2, 0 )
            vector<double> &new_out_parm,
            // whole_parm( 8, 0 )
            vector<double> whole_parm
            ) {
  int curr_label = ptr_seg[ 2 * ( idx - 1 ) ];
  if( sc != 0 ) {  // split or combine
    // energy of whole region, subregions and outlier
    int len = regions.size();
    vector<double> nrg( len + 1, 0 ); 
    vector<double>::iterator it_nrg = nrg.begin();
    // parameters for regions
    list<vector<double>> region_parm;
    // parameters for outler
    for( list<map<int, int>>::iterator it = regions.begin();
         it != regions.end(); ++ it, ++ it_nrg ) {
      double mu = - 1, sigma2 = 1;
      map<int, int>::iterator it_map = it->begin();
      int curr_label = it_map->second;
      updateParm( mu, theta, sigma2, *it, ptr_m[ 3 ], ptr_m[ 2 ],
                  ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst, curr_label,
                  ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, ptr_nintst,
                  ptr_alpha[ 3 ], ptr_beta[ 3 ], 20 );
      // curr_label, mu, sigma2, theta
      tmp_parm[ 0 ] = curr_label;
      tmp_parm[ 1 ] = mu;
      tmp_parm[ 2 ] = sigma2;
      for( int i = 0; i < 6; ++ i ) {
        tmp_parm[ i + 3 ] = theta[ i ];
      }

      region_parm.push_back( tmp_parm );

      // calculate energy
      set<int> t_region;
      for( ; it_map != it->end(); ++ it_map ) {
        t_region.insert( it_map->first );
      }
      double energy = energyY( t_region, mu, ptr_m[ 2 ], sigma2,
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
    map<int, int> out_region;
    out_region[ idx ] = out_label;
    updateParm( out_mu, out_sigma2, out_region, ptr_m[ 3 ], ptr_m[ 2 ],
                ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst, out_label,
                ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, ptr_alpha[ 3 ],
                ptr_beta[ 3 ], 20 );
    outlier_parm[ 0 ] = out_label;
    outlier_parm[ 1 ] = out_mu;
    outlier_parm[ 2 ] = out_sigma2;
    
    double out_energy = energyX( out_label, idx, true, ptr_seg, ptr_nidx,
                                 ptr_delta[ 0 ], ptr_gamma[ 0 ] );
    out_energy += energyY( out_label, idx, out_mu, ptr_m[ 2 ], out_sigma2,
                           ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, ptr_intst,
                           ptr_nintst, out_theta, ptr_alpha[ 3 ],
                           ptr_beta[ 3 ], ptr_a[ 0 ], ptr_b[ 0 ] );
    
    // single voxel energy ( -1, -2, -3 )
    double energy;
    double min_energy = out_energy;
    int min_label = out_label;
    
    double mu;
    double sigma2;
    for( int i = - 1; i > - 4; -- i ) {
      vector<double> &parm = health_parm[ i ];
      mu = parm[ 0 ];
      sigma2 = parm[ 1 ];
      for( int j = 0; j < 6; ++ j ) {
        theta[ j ] = parm[ 2 + j ];
      }
      energy = energyY( i, idx, mu, sigma2, ptr_seg, ptr_nidx,
                        ptr_intst, ptr_nintst, theta );
      energy += energyX( i, idx, true, ptr_seg, ptr_nidx,
                         ptr_delta[ 0 ], ptr_gamma[ 0 ] );
      if( energy < min_energy ) {
        min_energy = energy;
        min_label = i;
      }
    }
    
    *it_nrg = min_energy;
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
      // remove old whole region and add new splitted regions and
      // new outlier if min_label > 1
    } else if ( combine_nrg > split_nrg && sc == 1 ) {
      eraseRegion( whole_label, tumor_labels, tumor_regions, tumor_parm );
      list<vector<double>> new_region_parm( ++ region_parm.begin() ,
                                            region_parm.end() );
      list<map<int, int >> new_regions( ++ regions.begin(), 
                                        regions.end() );
      addRegion( ptr_seg, new_region_parm, new_regions, tumor_labels,
                 tumor_regions, tumor_parm );
      
      ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
      if( min_label > 0 ) {
        for( int i = 0; i < 2; ++ i ) {
          new_out_parm[ i ] = outlier_parm[ 1 + i ];
        }
        addOutl( out_label, new_out_parm, outl_labels, outl_parm );
      }
      // remove old subregions, parameters and labels,
      // add new whole regions label,region and parameters,
      // set seg to tumor,
      // remove outlier from outl_labels
      // remove outl_parm
    } else if( combine_nrg < split_nrg && sc == 2 ) {
      for( list<vector<double>>::iterator it_list =  ++ region_parm.begin();
           it_list != region_parm.end(); ++ it_list ) {
        int sub_label = it_list->front();
        eraseRegion( sub_label, tumor_labels, tumor_regions, tumor_parm );
      }
      
      list<vector<double>> new_whole_parm;
      list<map<int,int>> new_whole_region;
      new_whole_parm.push_back( label_whole_parm );
      new_whole_region.push_back( regions.front() );
      addRegion( ptr_seg, new_whole_parm, new_whole_region, tumor_labels,
                 tumor_regions, tumor_parm );
      ptr_seg[ 2 * ( idx - 1 ) ] = whole_label;
      if( curr_label > 0 ) {
        eraseOutl( out_label, outl_labels, outl_parm );
      }
      // update parameters for subregions
      // possibly remove or add outlier;
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
      ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
      // healthy to outlier
      if( min_label > 0 && curr_label < 1 ) {
        for( int i = 0; i < 2; ++ i ) {
          new_out_parm[ i ] = outlier_parm[ 1 + i ];
        }
        addOutl( out_label, new_out_parm, outl_labels, outl_parm );
        // outlier to healthy
      } else if( min_label < 0 && curr_label > 0 ) {
        eraseOutl( out_label, outl_labels, outl_parm );
      }
    }
    // no split or combine
  } else {
    // check if having tumor neighbor
    vector<int> nbr_label = nbrLabel( idx, ptr_seg, ptr_nidx );
    bool have_tumor = false;
    int t_label = 0;
    for( int i = 0; i < 6; ++ i ) {
      if( nbr_label[ i ] != NA_INTEGER && nbr_label[ i ] < - 3 ) {
        have_tumor = true;
        t_label = nbr_label[ i ];
      }
    }
    // have tumor neighbor
    if( have_tumor ) {
      // tumor energy
      vector<double> &t_parm = tumor_parm[ t_label ];
      double t_mu = t_parm[ 0 ];
      double t_sigma2 = t_parm[ 1 ];
      vector<double> &t_theta = theta;
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
      double min_energy = t_energy;
      int min_label = t_label;
      // Rprintf( "t_label = %d; t_energy = %f; t_mu = %f; t_sigma2 = %f;\n",
      //          t_label, t_energy, t_mu, t_sigma2 );
      // healthy cell
      double h_energy;
      
      double h_mu;
      double h_sigma2;
      vector<double> &h_theta = theta;
      for( int i = - 1; i > - 4; -- i ) {
        vector<double> &h_parm = health_parm[ i ];
        h_mu = h_parm[ 0 ];
        h_sigma2 = h_parm[ 1 ];
        for( int j = 0; j < 6; ++ j ) {
          h_theta[ j ] = h_parm[ 2 + j ];
        }
        h_energy = energyY( i, idx, h_mu, h_sigma2, ptr_seg, ptr_nidx,
                            ptr_intst, ptr_nintst, h_theta );
        h_energy += energyX( i, idx, false, ptr_seg, ptr_nidx,
                             ptr_delta[ 0 ], ptr_gamma[ 0 ] );
        // Rprintf( "h_label = %d; h_energy = %f; h_mu = %f; h_sigma2 = %f;\n",
        //          i, h_energy, h_mu, h_sigma2 );
        if( h_energy < min_energy ) {
          min_energy = h_energy;
          min_label = i;
        }
      }
      // outlier energy
      // outlier label
      int out_label = findOutLabel( idx, ptr_seg, outl_labels );
      
      
      double out_energy = energyX( out_label, idx, false, ptr_seg,
                                   ptr_nidx, ptr_delta[ 0 ],
                                   ptr_gamma[ 0 ] );
      // outlier parameters
      double out_mu = - 1, out_sigma2 = 1;
      map<int, int> out_region;
      out_region[ idx ] = out_label;
      if( curr_label > 0 ) {
        new_out_parm = outl_parm[ curr_label ];
        out_mu = new_out_parm[ 0 ];
        out_sigma2 = new_out_parm[ 1 ];
      } else {
        updateParm( out_mu, out_sigma2, out_region, ptr_m[ 3 ],
                    ptr_m[ 2 ], ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst,
                    out_label, ptr_lambda2[ 3 ], ptr_seg, ptr_nidx,
                    ptr_alpha[ 3 ], ptr_beta[ 3 ], 20 );
        new_out_parm[ 0 ]= out_mu;
        new_out_parm[ 1 ] = out_sigma2;
      }
      out_energy += energyY( out_label, idx, out_mu, ptr_m[ 2 ],
                             out_sigma2, ptr_lambda2[ 3 ], ptr_seg,
                             ptr_nidx, ptr_intst, ptr_nintst, out_theta,
                             ptr_alpha[ 3 ],ptr_beta[ 3 ], ptr_a[ 0 ],
                             ptr_b[ 0 ] );
      // Rprintf( "out_label = %d; out_energy = %f; out_mu = %f; out_sigma2 = %f\n",
      //          out_label, out_energy, out_mu, out_sigma2 );
      if( out_energy < min_energy ) {
        min_label = out_label;
      }
      if( min_label <= - 4 ) {
        if( curr_label >= - 3 && curr_label <= 0 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          tumor_regions[ t_label ].insert( idx );
        } else if( curr_label >= 1 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          eraseOutl( out_label, outl_labels, outl_parm );
          tumor_regions[ t_label ].insert( idx );
        }
      } else if( min_label >= - 3 && min_label <= -1 ) {
        if( curr_label == 0 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
        } else if( curr_label <= - 4 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          tumor_regions[ t_label ].erase( idx );
        } else if( curr_label >= 1 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          eraseOutl( out_label, outl_labels, outl_parm );
        }
      } else if( min_label >= 1 ) {
        if( curr_label  <= - 4 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          tumor_regions[ t_label ].erase( idx );
          addOutl( out_label, new_out_parm, outl_labels, outl_parm );
        } else if( curr_label >= - 3 && curr_label <= 0 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          addOutl( out_label, new_out_parm, outl_labels, outl_parm );
        }
      }
      
      // doesn't have tumor neighbor
    } else {
      // tumor energy
      
      if( curr_label <= - 4 ) {
        t_label = curr_label;
      } else {
        // Find a new tumor region label;
        t_label = newTumorLabel( 1, tumor_labels );
      }
      // new tumor region parameters
      double t_mu = -1, t_sigma2 = 1; // t_sigma2 has to be non-zero;
      vector<double> &t_theta = theta;
      map<int, int> new_region;
      list<map<int, int>> new_regions;
      new_region[ idx ] = t_label;
      new_regions.push_back( new_region );
      
      updateParm( t_mu, t_theta, t_sigma2, new_region, ptr_m[ 3 ],
                  ptr_m[ 2 ], ptr_a[ 0 ], ptr_b[ 0 ],
                  ptr_intst, t_label, ptr_lambda2[ 3 ], ptr_seg,
                  ptr_nidx, ptr_nintst, ptr_alpha[ 3 ],
                  ptr_beta[ 3 ], 20 );
      vector<double> &new_parm = tmp_parm;
      list<vector<double>> new_region_parm;
      new_parm[ 0 ] = t_label;
      new_parm[ 1 ] = t_mu;
      new_parm[ 2 ] = t_sigma2;
      for( int i = 0; i < 6; ++ i ) {
        new_parm[ i + 3 ] = t_theta[ i ];
      }
      new_region_parm.push_back( new_parm );
      
      double t_energy = energyY( t_label, idx, t_mu, ptr_m[ 2 ],
                                 t_sigma2, ptr_lambda2[ 3 ], ptr_seg,
                                 ptr_nidx, ptr_intst, ptr_nintst,
                                 t_theta, ptr_alpha[ 3 ], ptr_beta[ 3 ],
                                 ptr_a[ 0 ], ptr_b[ 0 ] );
      t_energy += energyX( t_label, idx, false, ptr_seg, ptr_nidx,
                           ptr_delta[ 0 ], ptr_gamma[ 0 ] );
      // Rprintf( "t_label = %d; t_energy = %f; t_mu = %f; t_sigma2 = %f\n", 
      //          t_label, t_energy, t_mu, t_sigma2 );
      double min_energy = t_energy;
      int min_label = t_label;
      
      // healthy cell
      double h_energy;
      
      double h_mu;
      double h_sigma2;
      vector<double> &h_theta = theta;
      for( int i = - 1; i > - 4; -- i ) {
        vector<double> &h_parm = health_parm[ i ];
        h_mu = h_parm[ 0 ];
        h_sigma2 = h_parm[ 1 ];
        for( int j = 0; j < 6; ++ j ) {
          h_theta[ j ] = h_parm[ 2 + j ];
        }
        h_energy = energyY( i, idx, h_mu, h_sigma2, ptr_seg, ptr_nidx,
                            ptr_intst, ptr_nintst, h_theta );
        h_energy += energyX( i, idx, false, ptr_seg, ptr_nidx,
                             ptr_delta[ 0 ], ptr_gamma[ 0 ] );
        // Rprintf( "h_label = %d; h_energy = %f; h_mu = %f; h_sigma2 = %f\n",
        //          i, h_energy, h_mu, h_sigma2 );
        if( h_energy < min_energy ) {
          min_energy = h_energy;
          min_label = i;
        }
      }
      if( min_label <= - 4 ) {
        if( curr_label >= - 3 && curr_label <= 0 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          addRegion( ptr_seg, new_region_parm, new_regions, tumor_labels,
                     tumor_regions, tumor_parm );
        } else if( curr_label >= 1 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          eraseOutl( curr_label, outl_labels, outl_parm );
          addRegion( ptr_seg, new_region_parm, new_regions, tumor_labels,
                     tumor_regions, tumor_parm );
        }
      } else if( min_label >= - 3 && min_label <= - 1 ) {
        if( curr_label == 0 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
        } else if( curr_label <= - 4 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          eraseRegion( t_label, tumor_labels, tumor_regions, tumor_parm );
        } else if( curr_label >= 1 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          eraseOutl( curr_label, outl_labels, outl_parm );
        }
      }
    }
  }
}
