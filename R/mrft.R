# Markov random field model for brain tumor segmentation
mrft <- function( patient, out = "SEG", infolder = "N4ITK433Z",
                  ## Always four numbers for delta
                  delta = 
                    list( t1ce = c( 0, 1, 9, 9 ),
                          flair = c( 0, 0, 8, 8 ),
                          t2 = c( 6, 0, NA_real_, NA_real_ ),
                          fthr = c( 0, 0, 8, 0 ) ),
                  gamma = list( t1ce = 0.8,
                                flair = 0.8,
                                t2 = 0.6,
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
                  min_enh = 2000L,
                  max_prop_enh_enc = .1,
                  min_tumor = 20000L,
                  spread = 4, 
                  min_prop_tumor_nbr = 0.6 ) {
  segment( patient, out, infolder, delta, gamma, alpha,
           beta, lambda2, a, nu2, maxit )
  post( patient, out, infolder, delta, gamma, alpha,
        beta, lambda2, a, nu2, maxit, 
        min_enh, max_prop_enh_enc, 
        min_tumor, spread, min_prop_tumor_nbr )
}