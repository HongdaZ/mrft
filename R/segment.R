# segment the MR images 
# beta and nu2 influenced by normalization

segment <- function( patient, delta = 5 ^ 2, gamma = 1, 
                     alpha = rep( 10, 4 ),
                     beta = rep( 1, 4 ),
                     lambda2 = rep( 1 / ( 36 * 4 ) , 4 ), 
                     a = 5,
                     nu2 = rep( .25, 3 ), 
                     maxit = 50L ) {
  data <- readImage( patient )
  t1ce_data <- splitT1ce( data$t1ce, data$flair )
  
}
