# segment the MR images 
# beta and nu2 influenced by normalization

segment <- function( patient, 
                     delta = list( t1ce = c( 4 ^ 2 / 2, 4 ^ 2 / 2 ),
                                   flair = c( 4 ^ 2 / 2, 4 ^ 2 / 2 ),
                                   t1 = c( 4 ^ 2 / 2, 4 ^ 2 / 2 ) ), 
                     gamma = 1, 
                     alpha = rep( 10, 4 ),
                     beta = rep( 1, 4 ),
                     lambda2 = rep( 1 , 4 ), 
                     a = 5,
                     nu2 = rep( .25, 3 ), 
                     maxit = list( t1ce = 10L, flair = 1L, t2 = 1L ),
                     t1ce_factor = 18L,
                     flair_factor = 18L ) {
  images <- readImage( patient )
  t1ce_data <- splitT1ce3( images$t1ce, images$flair )
  m <- t1ce_data$m
  t1ce_model <- initEst( t1ce_data$label, t1ce_data$t1ce )
  ## estimate parameters of t1ce or FLAIR images without tumor
  system.time( t1ce_seg <- est3( t1ce_model, delta$t1ce, gamma,
                                 alpha[ 1 : 3 ], beta[ 1 : 3 ], 
                                 lambda2[ 1 : 3 ], 
                                 m, nu2[ 1 : 3 ], 40L ) )
  
  ## update beta
  sigma2 <- t1ce_seg$parm[ 3, ]
  beta[ 1 : 3 ] <- ( alpha[ 1 : 3 ] + 1 ) * ( sigma2 )
  t1ce_data <- split4( t1ce_data$t1ce, t1ce_seg, t1ce_factor )
  ## update m
  m <- t1ce_data$m
  b <- getB( m, a )
  t1ce_model <- initEst( t1ce_data$label, t1ce_data$intst )
  sink( '/home/hzhang/Documents/t1ce_output.txt' )
  system.time( t1ce_seg <- pred4( t1ce_model, delta$t1ce, 
                                  gamma, alpha, beta, 
                                  lambda2, a, b, m, nu2, maxit$t1ce ) )
  
  sink()
  ## split flair images to CSF & necrosis, white matter and grey matter
  flair_data <- splitFlair3( images$flair, t1ce_seg )
  m <- flair_data$m
  flair_model <- initEst( flair_data$label, flair_data$flair )
  ## estimate parameters of t1ce or FLAIR images without tumor
  flair_seg <- est3( flair_model, delta$flair, gamma, 
                     alpha[ 1 : 3 ], beta[ 1 : 3 ], lambda2[ 1 : 3 ],
                     m, nu2[ 1 : 3 ], 40L )
  ## update beta
  sigma2 <- flair_seg$parm[ 3, ]
  beta[ 1 : 3 ] <- ( alpha[ 1 : 3 ] + 1 ) * ( sigma2 )
  flair_data <- split4( flair_data$flair, flair_seg, flair_factor )
  ## update m
  m <- flair_data$m
  b <- getB( m, a )
  flair_model <- initEst( flair_data$label, flair_data$intst )
  sink( '/home/hzhang/Documents/flair_output.txt' )
  system.time( flair_seg <- pred4( flair_model, delta$flair,
                                   gamma, alpha, beta,
                                   lambda2, a, b, m, nu2, 
                                   maxit$flair ) )
  sink()
  
}