# Markov random field model for brain tumor segmentation
mrft <- function( patient, out = "SEG", infolder = "N4ITK433Z",
                  ## Always four numbers for delta
                  delta = 
                    list( t1ce = c( -1, 0, 6, 2 ),
                          flair = c( -1, 0, NA_real_, 2 ),
                          t2 = c( 1, 0, NA_real_, 2 ),
                          fthr = c( 0, 0, 4, 8 ) ),
                  delta_factor = 
                    list( t1ce = 1 / 4,
                          flair = 1 / 4,
                          t2 = 1 / 9 ),
                  gamma = list( t1ce = 1,
                                flair = 1,
                                t2 = 1,
                                fthr = 1 ),
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
                  min_enh = 2000L,
                  min_enh_enc = 1000L,
                  max_prop_enh_enc = .1,
                  min_tumor = 20000L,
                  spread_add = 4,
                  spread_rm = 3 ) {
  segment( patient, out, infolder, delta, delta_factor,
           gamma, alpha,
           beta, lambda2, a, nu2, maxit )
  post( patient, out, infolder, delta, delta_factor,
        gamma, alpha,
        beta, lambda2, a, nu2, maxit, 
        min_enh, min_enh_enc, max_prop_enh_enc, 
        min_tumor, spread_add, spread_rm )
}