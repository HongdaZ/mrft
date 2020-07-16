#include <R.h>
#include <Rinternals.h>

#include <queue>
#include <stack>
#include <vector>
#include <set>

#include "search.h"
#include "findRegion.h"
#include "helper.h"
#include "initRegion.h"
#include "scTrn.h"
#include "scPred.h"
#include "updateMu.h"
#include "updateTS.h"
#include "updateSigma.h"
#include "energyY.h"
#include "energyX.h"
#include "nbrLabel.h"
#include "updateParm.h"
#include "initParm.h"
#include "cmpET.h"
#include "cmpEP.h"
#include "updateBeta.h"
#include "skip.h"
#include "debug.h"

using std::stack;
using std::queue;
using std::vector;
using std::set;

extern "C" {
  SEXP estParm( SEXP model, SEXP delta, SEXP gamma, 
                SEXP alpha, SEXP beta, SEXP lambda2, 
                SEXP a, SEXP b, SEXP m, SEXP nu2 ) {
    SEXP info = getListElement( model, "info" );
    SEXP seg = getListElement( model, "seg" );
    
    int *ptr_seg = INTEGER( seg );
    
    SEXP idx = getListElement( info, "idx" );
    SEXP nidx = getListElement( info, "nidx" );
    SEXP intst = getListElement( info, "intst" );
    SEXP nintst = getListElement( info, "nintst" );
    
    const int *ptr_idx = INTEGER( idx );
    const int *ptr_nidx = INTEGER( nidx );
    const double *ptr_intst = REAL( intst );
    const double *ptr_nintst = REAL( nintst );
    
    const double *ptr_delta = REAL( delta );
    const double *ptr_gamma = REAL( gamma );
    const double *ptr_alpha = REAL( alpha );
    double *ptr_beta = REAL( beta );
    const double *ptr_lambda2 = REAL( lambda2 );
    const double *ptr_a = REAL( a );
    const double *ptr_b = REAL( b );
    const double *ptr_m = REAL( m );
    const double *ptr_nu2 = REAL( nu2 );
    
    int len = length( idx );
    list<int> tumor_labels;
    list<int> outl_labels;
    map<int, vector<double>> health_parm;
    map<int, vector<double>> tumor_parm;
    map<int, vector<double>> outl_parm;

    
    map<int, list<int>> tumor_regions;
    initRegion( ptr_seg, ptr_nidx, len,
                tumor_regions, tumor_labels );
    
    // list< map<int, int>> regions;
    // // For testing scTrn
    // int flag = scTrn( regions, tumor_labels, tumor_regions, ptr_seg,
    //                   ptr_nidx, 1032015 );
    // Rprintf( "combine or split: %d \n", flag );
    // ///////////////////////////////////////////////////////////////
    // // Debug updateMu
    // for( int j = 1; j < 4; ++ j ) {
    //   set<int> region_h;
    //   for( int i = 0; i < len; ++ i ) {
    //     if( ptr_seg[ 2 * i ] == - j ) {
    //       region_h.insert( i + 1 ); // region starts from 1
    //     }
    //   }
    //   double sigma2 = 2;
    //   double m = 3;
    //   double nu2 = 2;
    //   vector<double> theta;
    //   for( int k = 0; k < 6; ++ k ) {
    //     theta.push_back( .01 + k * .01 );
    //   }
    //   double mu = updateMu( region_h, sigma2, m, nu2, theta, ptr_intst );
    //   Rprintf( " region %d, mu = %f", j, mu );
    // }
    //////////////////////////////////////////////////////////////////
    // // Debug updateMu for tumor
    // map<int, int> tumor_34;
    // for( int i = 0; i < len; ++ i ) {
    //   if( ptr_seg[ 2 * i ] == -34 ) {
    //     tumor_34[ i + 1 ] = -34;
    //   }
    // }
    // vector<double> theta;
    // for( int k = 0; k < 6; ++ k ) {
    //   theta.push_back( .01 + k * .01 );
    // }
    // double mu = updateMu( tumor_34, 2, 1, 2, 5, 6, theta, ptr_intst );
    // Rprintf( " Tumor region -34, mu = %f", mu );
    // // Debug updateTS
    // set<int> region_h;
    // for( int i = 0; i < len; ++ i ) {
    //   if( ptr_seg[ 2 * i ] == - 23 ) {
    //     region_h.insert( i + 1 ); // region starts from 1
    //   }
    // }
    // double sigma2 = 2;
    // vector<double> theta;
    // updateTS( region_h, -23, 1, sigma2, 3, ptr_seg, ptr_nidx, ptr_intst, 
    //           ptr_nintst, theta, 1, .00001 );
    // // updateSigma( 20, 1, ptr_intst, 2, .5, sigma2 );
    // // Rprintf( "sigma2 = %f\n", sigma2 );
    // // for( int k = 0; k < 6; ++ k ) {
    // //   Rprintf( "theta = %f\t", theta[ k ] );
    // // }
    // // Rprintf( "\n" );
    // ////////////////////////////////////////////////////////////////////////
    // double energy = energyY( region_h, 1, .5, sigma2, 3, ptr_seg, ptr_nidx,
    //                          ptr_intst, ptr_nintst, theta, 
    //                          1, .00001, ptr_a[ 0 ], ptr_b[ 0 ] );
    // // Rprintf( "energyY = %f\n", energy );
    // vector<double> theta;
    // double vtheta[ 6 ] = { 0.1, 0.3, 0.4, 0.2, 0.5, 6.0 };
    // for( int i = 0; i < 6; ++ i ) {
    //   theta.push_back( vtheta[ i ] ) ;
    // }
    // // double energy = energyY( -1, 1208, 1, 0.13, ptr_seg, 
    // //                          ptr_nidx, ptr_intst, ptr_nintst, theta );
    // // Rprintf( "energyY = %f\n", energy);
    // double energy = energyY( 1, 504138, 1, .5, 0.13, 3, ptr_seg, 
    //                          ptr_nidx, ptr_intst, ptr_nintst, theta,
    //                          1, .00001, ptr_a[ 0 ], ptr_b[ 0 ] );
    // // Rprintf( "energyY = %f\n", energy);
    // 
    // 
    // energy = energyX( 1, 507521, true,
    //                 ptr_seg, ptr_nidx, ptr_delta[ 0 ], ptr_gamma[ 0 ] );
    // // Rprintf( "energyX = %f\n", energy);
    // vector<int> nbr_label = nbrLabel( 1059, ptr_seg, ptr_nidx );
    // for( int i = 0; i < 6; ++ i ) {
    //   Rprintf( "neighbor label:%d\n", nbr_label[ i ] );
    // }
    // // debug updateParm
    // set<int> region_healthy;
    // int curr_label = - 2;
    // for( int i = 0; i < len; ++ i ) {
    //   if( ptr_seg[ 2 * i ] == curr_label ) {
    //     region_healthy.insert( i + 1 ); // region starts from 1
    //   }
    // }
    // double mu = 1.0;
    // double sigma2 = 2.0;
    // updateParm( mu, theta, sigma2, region_healthy, ptr_m[ 1 ], ptr_nu2[ 1 ],
    //             ptr_intst, curr_label, ptr_lambda2[ 1 ], ptr_seg, ptr_nidx,
    //             ptr_nintst, ptr_alpha[ 1 ], ptr_beta[ 1 ], 20 );
    ////////////////////////////////////////////////////////////////////
    // // debug updateParm
    // map<int, int> region_t;
    // int curr_label = - 34;
    // for( int i = 0; i < len; ++ i ) {
    //   if( ptr_seg[ 2 * i ] == curr_label ) {
    //     region_t[ i + 1 ] = curr_label; // region starts from 1
    //   }
    // }
    // double mu = 1.0;
    // double sigma2 = 2.0;
    // updateParm( mu, theta, sigma2, region_t, ptr_m[ 3 ], ptr_m[ 2 ], 
    //             ptr_a[ 0 ],
    //             ptr_b[ 0 ], ptr_intst, curr_label, ptr_lambda2[ 3 ], 
    //             ptr_seg,
    //             ptr_nidx, ptr_nintst, ptr_alpha[ 3 ], ptr_beta[ 3 ], 20 );
    /////////////////////////////////////////////////////////////////////////
    // // Debug updateParm for outliers
    // int curr_label = 1;
    // map<int, int> region_o;
    // region_o[ 1060076 ] = curr_label;
    // double mu = 1.0;
    // double sigma2 = 2.0;
    // updateParm( mu, sigma2, region_o, ptr_m[ 3 ], ptr_m[ 2 ], ptr_a[ 0 ], 
    //             ptr_b[ 0 ], ptr_intst, curr_label, ptr_lambda2[ 3 ], 
    //             ptr_seg, ptr_nidx, ptr_alpha[ 3 ], ptr_beta[ 3 ], 20 );
    
    initParm( health_parm, tumor_parm, ptr_seg, ptr_m, ptr_nu2, ptr_intst,
              ptr_lambda2, ptr_nidx, ptr_nintst, ptr_alpha, ptr_beta,
              tumor_regions, ptr_a, ptr_b, len, 20 );
    // // update beta
    // for( int i = 0; i < 3; ++ i ){
    //   double sigma2 = health_parm[ - i - 1 ][ 1 ];
    //   ptr_beta[ i ] = ( ptr_alpha[ i ] + 1 ) * sigma2;
    // }
    // // update beta for tumor regions based on the date of 
    // // the biggest tumor region
    // int max_size = 0;
    // int max_label = 0;
    // for( map<int, set<int>>::iterator it = tumor_regions.begin();
    //      it != tumor_regions.end(); ++ it ) {
    //   if( it->second.size() > max_size ) {
    //     max_size = it->second.size();
    //     max_label = it->first;
    //   }
    // }
    // double tumor_sigma2 = tumor_parm[ max_label ][ 1 ];
    // ptr_beta[ 3 ] = ( ptr_alpha[ 3 ] + 1 ) * tumor_sigma2;
    updateBeta( ptr_beta, ptr_alpha, health_parm, tumor_regions, 
                tumor_parm );
    // for( int i = 0; i < 4; ++ i ) {
    //   Rprintf( "beta = %f\t", ptr_beta[ i ] );
    // }
    // Rprintf( "\n" );
    // for( map<int, vector<double>>::iterator it = health_parm.begin();
    //      it!= health_parm.end(); ++ it) {
    //   Rprintf( "label = %d:", it->first );
    //   for( int i = 0; i < 8; ++ i ) {
    //     Rprintf( "\t%f", it->second[ i ] );
    //   }
    //   Rprintf( "\n" );
    // }
    // for( map<int, vector<double>>::iterator it = tumor_parm.begin(); 
    //      it!= tumor_parm.end(); ++ it) {
    //   Rprintf( "label = %d", it->first );
    //   for( int i = 0; i < 8; ++ i ) {
    //     Rprintf( "\t%f", it->second[ i ] );
    //   }
    //   Rprintf( "\n" );
    // }
    // 
    // 
    // // debug updateParm
    // set<int> region_healthy;
    // int curr_label = - 2;
    // for( int i = 0; i < len; ++ i ) {
    //   if( ptr_seg[ 2 * i ] == curr_label ) {
    //     region_healthy.insert( i + 1 ); // region starts from 1
    //   }
    // }
    // double mu = 1.0;
    // double sigma2 = 2.0;
    // vector<double> theta;
    // for( int i = 0; i < 6; ++ i ) {
    //   theta.push_back( 0 );
    // }
    // updateParm( mu, theta, sigma2, region_healthy, ptr_m[ 1 ], ptr_nu2[ 1 ],
    //             ptr_intst, curr_label, ptr_lambda2[ 1 ], ptr_seg, ptr_nidx,
    //             ptr_nintst, ptr_alpha[ 1 ], ptr_beta[ 1 ], 20 );
    // 
    // Rprintf( "label = %d\t%f\t%f\t", curr_label, mu, sigma2 );
    // for( int i = 0; i < 6; ++ i ) {
    //   Rprintf( "%f\t", theta[ i ] );
    // }
    // ///////////////////////////////////////////////////////////////////////
    
    // tumor_labels.erase( -42 );
    // list< map<int, int>> regions; //sub-regions
    // int flag = scTrn( regions, tumor_labels, tumor_regions, ptr_seg,
    //                   ptr_nidx, 1032015 );
    // Rprintf( "size = %d\n", regions.size() );
    // for( list<map<int, int>>::iterator it = regions.begin();
    //      it != regions.end(); ++ it ) {
    //   map<int, int> region = *it;
    //   for( map<int, int>::iterator it_map = region.begin(); 
    //        it_map != region.end(); ++ it_map ) {
    //     Rprintf( "%d = %d;", it_map->first, it_map->second );
    //   }
    //   Rprintf( "\n" );
    // }
    // list< map<int, int>> regions;
    // int curr_idx = 1032015;1032015
    // int curr_idx = 1032012;
    // int sc = scTrn( regions, tumor_labels, tumor_regions, ptr_seg,
    //                   ptr_nidx, curr_idx );
    // cmpET( curr_idx, sc, regions, tumor_regions, tumor_labels, outl_labels,
    //        health_parm, tumor_parm, outl_parm, ptr_seg, ptr_nidx, ptr_intst,
    //        ptr_nintst, ptr_delta, ptr_gamma, ptr_alpha, ptr_beta,
    //        ptr_lambda2, ptr_a, ptr_b, ptr_m, ptr_nu2 );
    // // debug scPred
    // curr_idx = 1032015;
    // sc = scPred( regions, tumor_labels, tumor_regions, ptr_seg,
    //              ptr_nidx, curr_idx );
    // Rprintf( "sc = %d\n", sc );
    // for( list<map<int, int>>::iterator it_list = regions.begin();
    //      it_list != regions.end(); ++ it_list ) {
    //   Rprintf( "\n" );
    //   for( map<int, int>::iterator it_map = it_list->begin();
    //        it_map != it_list->end(); ++ it_map ) {
    //     Rprintf( "idx = %d; label = %d; ", it_map->first, it_map->second );
    //   }
    // }
    // //////////////
    
    // list< map<int, int>> regions;
    // int curr_idx = 542755;
    // int sc = scPred( regions, tumor_labels, tumor_regions, ptr_seg,
    //                   ptr_nidx, curr_idx );
    // // Rprintf( "sc = %d\n", sc );
    // vector<double> outlier_parm( 3, 0 ), theta( 6, 0 ), tmp_parm( 9, 0 );
    // vector<double> out_theta( 6, 0 ), new_out_parm( 2, 0 );
    // vector<double> whole_parm( 8, 0 );
    // cmpEP( curr_idx, sc, regions, tumor_regions, tumor_labels, outl_labels,
    //        health_parm, tumor_parm, outl_parm, ptr_seg, ptr_nidx, ptr_intst,
    //        ptr_nintst, ptr_delta, ptr_gamma, ptr_alpha, ptr_beta,
    //        ptr_lambda2, ptr_a, ptr_b, ptr_m, ptr_nu2, outlier_parm,
    //        theta, tmp_parm, out_theta, new_out_parm, whole_parm );
    // return seg;
    
    //debug();
    
    int maxit = 20;
    bool skip_curr;
    vector<int> search( len, 0 );
    double threshold = ptr_m[ 2 ];
    for( int i = 0; i < maxit; ++ i ) {
      for( int j = 1; j <= len; ++ j ) {
        
        skip_curr = skip( j, ptr_seg, ptr_nidx, ptr_intst, threshold, 
                          search[ j ], 10 );
        if( skip_curr ) {
          continue;
        } else {
          ++ search[ j ];
          list<list<int>> regions;
          list<int> labels;
          int sc = scTrn( labels, regions, tumor_labels, tumor_regions, 
                          ptr_seg, ptr_nidx, j );
          cmpET( j, sc, labels, regions, tumor_regions, tumor_labels, 
                 outl_labels, health_parm, tumor_parm, outl_parm, ptr_seg,
                 ptr_nidx, ptr_intst, ptr_nintst, ptr_delta, ptr_gamma,
                 ptr_alpha, ptr_beta, ptr_lambda2, ptr_a, ptr_b, ptr_m,
                 ptr_nu2 );
        }
        // Rprintf( "%d\t; curr_label = %d\n", j, ptr_seg[ 2 * ( j - 1 ) ] );
      }
      // update parm for healthy regions
      for( int i = - 1; i > - 4;  -- i ) {
        int curr_label = i;
        // Rprintf( "curr_label = %d \n", curr_label );
        double mu = -1, sigma2 = 1;
        vector<double> theta;
        for( int  j = 0; j < 6; ++ j ) {
          theta.push_back( 0 );
        }
        list<int> region;
        for( int k = 0; k < len; ++ k ) {
          if( ptr_seg[ 2 * k ] == curr_label ) {
            region.push_back( k + 1 ); // region starts from 1
          }
        }
        int h_idx =  - 1 - curr_label; // == 0, 1, 2
        updateParm( mu, theta, sigma2, region, ptr_m[ h_idx ], ptr_nu2[ h_idx ],
                    ptr_intst, curr_label, ptr_lambda2[ h_idx ], ptr_seg, ptr_nidx,
                    ptr_nintst, ptr_alpha[ h_idx ], ptr_beta[ h_idx ], maxit );

        vector<double> h_parm;
        h_parm.push_back( mu );
        h_parm.push_back( sigma2 );

        h_parm.insert( h_parm.end(), theta.begin(), theta.end() );
        health_parm[ curr_label ] = h_parm;
      }
      // update parm for tumor_regions
    }
    return seg;
  }
} // extern "C"
