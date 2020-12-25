# Markov random field model for brain tumor segmentation
mrft <- function( patient, out = "SEG", infolder = "N4ITK433Z",
                  ## Always four numbers for delta
                  delta = 
                    list( t1ce = c( -1, 0, 5, 4 ),
                          flair = c( -0.5, 0, NA_real_, 4 ),
                          t2 = c( 0.5, 0, NA_real_, 4 ),
                          fthr = c( 0, 0, 4, 5 ) ),
                  delta_factor = 
                    list( t1ce = 1.75,
                          flair = 2.70,
                          t2 = 4.50 ),
                  gamma = list( t1ce = 0.8,
                                flair = 0.5,
                                t2 = 0.8,
                                fthr = 0.8 ),
                  ## #of healthy tissue types controlled by alpha
                  alpha = list( t1ce = rep( 10, 4 ),
                                flair = rep( 10, 4 ),
                                t2 = rep( 10, 3 ),
                                fthr = rep( 10, 2 ) ),
                  beta = list( t1ce = rep( 1, 4 ),
                               flair = rep( 1, 4 ),
                               t2 = rep( 1, 3 ),
                               fthr = rep( 1, 2 ) ),
                  lambda2 = list( t1ce = rep( 1 , 4 ),
                                  flair = rep( 1 , 4 ),
                                  t2 = rep( 1 , 3 ),
                                  fthr = rep( 1, 2 ) ),
                  a = 5,
                  nu2 = list( t1ce = rep( .25, 3 ),
                              flair = rep( .25, 3 ),
                              t2 = c( .25, .25 ),
                              fthr = rep( .25, 2 ) ),
                  maxit = 
                    list( t1ce = 10L, flair = 1L,
                          t2 = 1L, fthr = 40L ),
                  min_enh = 500L,
                  min_enh_enc = 2000L,
                  max_prop_enh_enc = .1,
                  max_prop_enh_slice = .2,
                  min_tumor = 20000L,
                  spread_add = 3.5,
                  spread_rm = 3.5,
                  spread_trim = 4 ) {
  infile <- patient[ 1 ]
  outfile <- gsub( infolder, out, infile )
  out_new_delta_t2 <- gsub( "_flair.nii.gz", "_post.rds", outfile )
  redo <- TRUE
  if( file.exists( out_new_delta_t2 ) ) {
    new_delta_t2 <- readRDS( out_new_delta_t2 )
    if( is.null( new_delta_t2 ) ) {
      redo <- FALSE
    } else {
      if( new_delta_t2[ 3 ] < 2 ) {
        redo <- FALSE
      } else {
        delta$t2 <- new_delta_t2
      }
    }
  }
  while( redo ) {
    segment( patient, out, infolder, delta, delta_factor,
             gamma, alpha,
             beta, lambda2, a, nu2, maxit, redo )
    post( patient, out, infolder, delta, delta_factor,
          gamma, alpha,
          beta, lambda2, a, nu2, maxit, 
          min_enh, min_enh_enc, max_prop_enh_enc, 
          max_prop_enh_slice,
          min_tumor, spread_add, spread_rm, spread_trim )
    if( file.exists( out_new_delta_t2 ) ) {
      new_delta_t2 <- readRDS( out_new_delta_t2 )
      if( is.null( new_delta_t2 ) ) {
        redo <- FALSE
      } else {
        if( new_delta_t2[ 3 ] < 2 ) {
          redo <- FALSE
        } else {
          delta$t2 <- new_delta_t2
        }
      }
    }
  }
}
