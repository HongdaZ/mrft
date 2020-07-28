# segment the MR images 
# beta and nu2 influenced by normalization

segment <- function( patient, delta = 5 ^ 2, gamma = 1, 
                     alpha = rep( 10, 4 ),
                     beta = rep( 1, 4 ),
                     lambda2 = rep( 1 / ( 36 * 4 ) , 4 ), 
                     a = 5,
                     nu2 = rep( .25, 3 ), 
                     maxit = 50L ) {
  images <- readImage( patient )
  t1ce_data <- splitT1ce3( images$t1ce, images$flair )
  m <- t1ce_data$m
  t1ce_model <- initEst( t1ce_data$label, t1ce_data$t1ce )
  # estimate parameters of t1ce or t2 images without tumor
  t1ce_seg <- est3( t1ce_model, delta, gamma,
                     alpha[ 1 : 3 ], beta[ 1 : 3 ], lambda2[ 1 : 3 ], 
                     m, nu2[ 1 : 3 ], 50L )
  t1ce_data <- splitT1ce4( images$t1ce, t1ce_seg )
  
  
  
}
