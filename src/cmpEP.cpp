#include <R.h>
#include <Rinternals.h>

#include "cmpEP.h"
#include "findOutLabel.h"
#include "newTumorLabel.h"
#include "eraseRegion.h"
#include "eraseOutl.h"
#include "addRegion.h"
#include "addOutl.h"
#include "tumorNbr.h"
#include "getRegion.h"
#include "assignParm.h"
#include "label2col.h"
#include "getParm.h"
#include "zeroVector.h"

// compare energy for prediction
void cmpEP( vector<int> &region, int idx, int sc,
            const vector<int> &labels, const vector<int> &regions,
            vector<int> &tumor_labels, vector<int> &outl_labels,
            const vector<double> &health_parm,
            vector<double> &tumor_parm, 
            vector<double> &outl_parm,
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
            int &n_tumor, int &n_outl,
            const int &len ) {
  int curr_label = ptr_seg[ 2 * ( idx - 1 ) ];
  if( sc != 0 ) {  // split or combine
    // energy of whole region, subregions and outlier or healthy cell
    int n_region = labels.size();
    vector<double> nrg( n_region + 1, 0 ); 
    // parameters for regions
    int nrow = 9;
    vector<double> region_parm( nrow * n_region, 0 );
    double mu;
    double sigma2;
    for( int i = 0; i < n_region; ++ i ) {
      int curr_label = labels[ i ];
      int row = ( i == 0 ) ? 0 : 1;
      getRegion( region, curr_label, regions, len, row );
      updateParm( mu, theta, sigma2, region, ptr_m[ 3 ], ptr_m[ 2 ],
                  ptr_a[ 0 ], ptr_b[ 0 ], ptr_intst, curr_label,
                  ptr_lambda2[ 3 ], ptr_seg, ptr_nidx, ptr_nintst,
                  ptr_alpha[ 3 ], ptr_beta[ 3 ], 20 );
      // curr_label, mu, sigma2, theta
      region_parm[ nrow * i ] = curr_label;
      region_parm[ nrow * i + 1 ] = mu;
      region_parm[ nrow * i + 2 ] = sigma2;
      for( int j = 0; j < 6; ++ j ) {
        region_parm[ nrow * i + j + 3 ] = theta[ j ];
      }
      // calculate energy
      double energy = energyY( region, mu, ptr_m[ 2 ], sigma2,
                               ptr_lambda2[ 3 ], ptr_seg, ptr_nidx,
                               ptr_intst, ptr_nintst,
                               theta, ptr_alpha[ 3 ], ptr_beta[ 3 ],
                               ptr_a[ 0 ], ptr_b[ 0 ] );
      nrg[ i ] = energy;
    }
    // New outlier label
    int out_label = findOutLabel( idx, ptr_seg, outl_labels );

    // outlier parameters
    double &out_mu = mu;
    double &out_sigma2 = sigma2;
    vector<int> out_region;
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
    out_energy += energyY( out_label, idx, out_mu, ptr_m[ 2 ], out_sigma2,
                           ptr_lambda2[ 3 ], ptr_seg, ptr_intst,
                           ptr_alpha[ 3 ],
                           ptr_beta[ 3 ], ptr_a[ 0 ], ptr_b[ 0 ] );
    
    // single voxel energy ( -1, -2, -3 )
    double energy;
    double min_energy = out_energy;
    int min_label = out_label;
    
    for( int i = - 1; i > - 4; -- i ) {
      int cidx = label2col( i );
      getParm( mu, sigma2, theta, health_parm, cidx );
      energy = energyY( i, idx, mu, sigma2, ptr_seg, ptr_nidx,
                        ptr_intst, ptr_nintst, theta );
      energy += energyX( i, idx, true, ptr_seg, ptr_nidx,
                         ptr_delta[ 0 ], ptr_gamma[ 0 ] );
      if( energy < min_energy ) {
        min_energy = energy;
        min_label = i;
      }
    }
    
    nrg.back() = min_energy;
    double combine_nrg = nrg[ 0 ];
    double split_nrg = 0;
    for( int i = 1; i < ( n_region + 1 ); ++ i ) {
      split_nrg += nrg[ i ];
    }
    vector<double> label_whole_parm( region_parm.begin(), 
                                     region_parm.begin() + nrow );
    int whole_label = label_whole_parm[ 0 ];
    // update parameters for whole region
    if( combine_nrg <= split_nrg && sc == 1 ) {
      for( int i = 0; i < 8; ++ i ) {
        whole_parm[ i ] = label_whole_parm[ 1 + i ];
      }
      assignParm( tumor_parm, whole_label, whole_parm );
      // remove old whole region and add new splitted regions and
      // new outlier if min_label > 1
    } else if ( combine_nrg > split_nrg && sc == 1 ) {
      eraseRegion( whole_label, tumor_labels, tumor_parm,
                   n_tumor );
      vector<double> new_region_parm( region_parm.begin() + nrow,
                                      region_parm.end() );
      addRegion( ptr_seg, new_region_parm, regions, 1, tumor_labels,
                 tumor_parm, n_tumor );
      
      ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
      if( min_label > 0 ) {
        for( int i = 0; i < 2; ++ i ) {
          new_out_parm[ i ] = outlier_parm[ 1 + i ];
        }
        addOutl( out_label, new_out_parm, outl_labels, outl_parm, n_outl );
      }
      // remove old subregions, parameters and labels,
      // add new whole regions label,region and parameters,
      // set seg to tumor,
      // remove outlier from outl_labels
      // remove outl_parm
    } else if( combine_nrg < split_nrg && sc == 2 ) {
      for( int i = 1; i < n_region; ++ i ) {
        int sub_label = labels[ i ];
        eraseRegion( sub_label, tumor_labels, tumor_parm, n_tumor );
      }
      addRegion( ptr_seg, label_whole_parm, regions, 0, tumor_labels,
                 tumor_parm, n_tumor );
      ptr_seg[ 2 * ( idx - 1 ) ] = whole_label;
      if( curr_label > 0 ) {
        eraseOutl( out_label, outl_labels, outl_parm, n_outl );
      }
      // update parameters for subregions
      // possibly remove or add outlier;
    } else if( combine_nrg >= split_nrg && sc == 2 ) {
      for( int i = 1; i < n_region; ++ i ) {
        int sub_label = labels[ i ];
        vector<double> &new_parm = whole_parm; 
        for( int j = 0; j < 8; ++ j ) {
          new_parm[ j ] = region_parm[ nrow * i + j + 1 ];
        }
        assignParm( tumor_parm, sub_label, new_parm );
      }
      ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
      // healthy to outlier
      if( min_label > 0 && curr_label < 1 ) {
        for( int i = 0; i < 2; ++ i ) {
          new_out_parm[ i ] = outlier_parm[ 1 + i ];
        }
        addOutl( out_label, new_out_parm, outl_labels, outl_parm, n_outl );
        // outlier to healthy
      } else if( min_label < 0 && curr_label > 0 ) {
        eraseOutl( out_label, outl_labels, outl_parm, n_outl );
      }
    }
    // no split or combine
  } else {
    // check if having tumor neighbor
    int tumor_idx = firstTumorNbr( idx, ptr_seg, ptr_nidx ) ;
    int t_label;
    double mu, sigma2;
    vector<double> &t_theta = theta;
    int cidx;
    // have tumor neighbor
    if( tumor_idx != 0 ) {
      int t_idx = tumor_idx;
      t_label = ptr_seg[ 2 * ( t_idx - 1 ) ];
      // tumor energy
      int cidx = label2col( t_label );
      double &t_mu = mu, &t_sigma2 = sigma2;
      getParm( t_mu, t_sigma2, t_theta, tumor_parm, cidx );
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
      
      double &h_mu = mu, &h_sigma2 = sigma2;
      vector<double> &h_theta = theta;
      for( int i = - 1; i > - 4; -- i ) {
        cidx = label2col( i );
        getParm( h_mu, h_sigma2, h_theta, health_parm, cidx );
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
      double &out_mu = mu, &out_sigma2 = sigma2;
      vector<int> out_region;
      out_region.push_back( idx );
      if( curr_label > 0 ) {
        cidx = label2col( curr_label );
        getParm( out_mu, out_sigma2, outl_parm, cidx );
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
                             ptr_intst,
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
        } else if( curr_label >= 1 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          eraseOutl( out_label, outl_labels, outl_parm, n_outl );
        }
      } else if( min_label >= - 3 && min_label <= -1 ) {
        if( curr_label == 0 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
        } else if( curr_label <= - 4 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
        } else if( curr_label >= 1 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          eraseOutl( out_label, outl_labels, outl_parm, n_outl );
        }
      } else if( min_label >= 1 ) {
        if( curr_label  <= - 4 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          addOutl( out_label, new_out_parm, outl_labels, outl_parm,
                   n_outl );
        } else if( curr_label >= - 3 && curr_label <= 0 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          addOutl( out_label, new_out_parm, outl_labels, outl_parm,
                   n_outl );
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
      // t_sigma2 has to be non-zero;
      double &t_mu = mu, &t_sigma2 = sigma2; 
      // use the function for outliers and single voxel tumor regions
      updateParm( t_mu, t_sigma2, idx, ptr_m[ 3 ],
                  ptr_m[ 2 ], ptr_a[ 0 ], ptr_b[ 0 ],
                  ptr_intst, ptr_seg,
                  ptr_alpha[ 3 ],
                  ptr_beta[ 3 ], 20 );
      vector<double> &new_parm = tmp_parm;
      new_parm[ 0 ] = t_label;
      new_parm[ 1 ] = t_mu;
      new_parm[ 2 ] = t_sigma2;
      for( int i = 0; i < 6; ++ i ) {
        new_parm[ i + 3 ] = t_theta[ i ];
      }
      double t_energy = energyY( t_label, idx, t_mu, ptr_m[ 2 ],
                                 t_sigma2, ptr_lambda2[ 3 ], ptr_seg,
                                 ptr_intst, ptr_alpha[ 3 ], ptr_beta[ 3 ],
                                 ptr_a[ 0 ], ptr_b[ 0 ] );
      t_energy += energyX( t_label, idx, false, ptr_seg, ptr_nidx,
                           ptr_delta[ 0 ], ptr_gamma[ 0 ] );
      // Rprintf( "t_label = %d; t_energy = %f; t_mu = %f; t_sigma2 = %f\n",
      //          t_label, t_energy, t_mu, t_sigma2 );
      double min_energy = t_energy;
      int min_label = t_label;
      
      // healthy cell
      double h_energy;
      double &h_mu = mu, &h_sigma2 = sigma2;
      vector<double> &h_theta = theta;
      for( int i = - 1; i > - 4; -- i ) {
        cidx = label2col( i );
        getParm( h_mu, h_sigma2, h_theta, health_parm, cidx );
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
          addRegion( ptr_seg, new_parm, new_region, tumor_labels,
                     tumor_parm, n_tumor );
        } else if( curr_label >= 1 ) {
          eraseOutl( curr_label, outl_labels, outl_parm, n_outl );
          addRegion( ptr_seg, new_parm, new_region, tumor_labels,
                     tumor_parm, n_tumor );
        }
      } else if( min_label >= - 3 && min_label <= - 1 ) {
        if( curr_label == 0 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
        } else if( curr_label <= - 4 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          eraseRegion( t_label, tumor_labels, tumor_parm,
                       n_tumor );
        } else if( curr_label >= 1 ) {
          ptr_seg[ 2 * ( idx - 1 ) ] = min_label;
          eraseOutl( curr_label, outl_labels, outl_parm, n_outl );
        }
      }
    }
  }
  return;
}
