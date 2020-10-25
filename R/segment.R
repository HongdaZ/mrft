# segment the MR images 
# beta and nu2 influenced by normalization

segment <- function( patient, 
                     delta = 
                       list( t1ce = c( 0, 0, 4 ^ 2 / 2, 4 ^ 2 / 2 ),
                             flair = c( 0, 0, 10, 4 ^ 2 / 2 ),
                             t2 = c( 8, 0, 4 ^ 2 / 2, 4 ^ 2 / 2 ),
                             fthr = c( 0, 0, 0, 0 ) ), 
                     gamma = 1, 
                     alpha = list( t1ce = rep( 10, 4 ),
                                   flair = rep( 10, 4 ),
                                   t2 = rep( 10, 3 ) ),
                     beta = list( t1ce = rep( 1, 4 ),
                                  flair = rep( 1, 4 ),
                                  t2 = rep( 1, 3 ) ),
                     lambda2 = list( t1ce = rep( 1 , 4 ),
                                     flair = rep( 1 , 4 ),
                                     t2 = rep( 1 , 3 ) ),
                     a = 5,
                     nu2 = list( t1ce = rep( .25, 3 ),
                                 flair = rep( .25, 3 ),
                                 t2 = c( .25, .01 ) ),
                     maxit = list( t1ce = 10L, 
                                   flair = 1L, t2 = 1L ),
                     factor = list( t1ce = 18L, 
                                    flair = 18L, t2 = 18L ) ) {
  images <- readImage( patient )
  ## split t1ce images to CSF & necrosis, grey matter and white matter
  t1ce_data <- splitT1ce3( images$t1ce, images$flair )
  m <- t1ce_data$m
  t1ce_model <- initEst( t1ce_data$label, t1ce_data$t1ce )
  ## estimate parameters of t2 images without tumor and CSF & necrosis
  t1ce_seg <- est( t1ce_model, delta$t1ce, gamma,
                    alpha$t1ce[ 1 : 3 ], beta$t1ce[ 1 : 3 ], 
                    lambda2$t1ce[ 1 : 3 ], 
                    m, nu2$t1ce, 40L )
  ## update beta
  sigma2 <- t1ce_seg$parm[ 3, ]
  beta$t1ce[ 1 : 3 ] <- ( alpha$t1ce[ 1 : 3 ] + 1 ) * ( sigma2 )
  t1ce_data <- split4( t1ce_data$t1ce, t1ce_seg, factor$t1ce )
  ## update m
  m <- t1ce_data$m
  b <- getB( m, a )
  t1ce_model <- initEst( t1ce_data$label, t1ce_data$intst )
  # sink( '/home/hzhang/Documents/t1ce_output.txt' )
  # system.time( t1ce_seg <- pred( t1ce_model, delta$t1ce, 
  #                                 gamma, alpha$t1ce, beta$t1ce, 
  #                                 lambda2$t1ce, a, b, m, nu2$t1ce,
  #                                 maxit$t1ce ) )
  # 
  # sink()
  t1ce_seg <- pred( t1ce_model, delta$t1ce, 
                     gamma, alpha$t1ce, beta$t1ce, 
                     lambda2$t1ce, a, b, m, nu2$t1ce, maxit$t1ce )
  ## Both tumor and outliers are regarded as tumor
  t1ce_seg$image[ t1ce_seg$image < -3 | t1ce_seg$image > 0 ] <- -4L
  
  ## split flair images to CSF & necrosis, white matter 
  ## and grey matter
  flair_data <- splitFlair3( images$flair, t1ce_seg )
  m <- flair_data$m
  flair_model <- initEst( flair_data$label, flair_data$flair )
  ## estimate parameters of t1ce or FLAIR images without tumor
  flair_seg <- est( flair_model, delta$flair, gamma, 
                     alpha$flair[ 1 : 3 ], beta$flair[ 1 : 3 ], 
                     lambda2$flair[ 1 : 3 ],
                     m, nu2$flair, 40L )
  ## update beta
  sigma2 <- flair_seg$parm[ 3, ]
  beta$flair[ 1 : 3 ] <- ( alpha$flair[ 1 : 3 ] + 1 ) * ( sigma2 )
  flair_data <- split4( flair_data$flair, flair_seg, factor$flair )
  ## update m
  m <- flair_data$m
  b <- getB( m, a )
  flair_model <- initEst( flair_data$label, flair_data$intst )
  # sink( '/home/hzhang/Documents/flair_output.txt' )
  # system.time( flair_seg <- pred( flair_model, delta$flair,
  #                                  gamma, alpha$flair, beta$flair,
  #                                  lambda2$flair, a, b, m, nu2$flair, 
  #                                  maxit$flair ) )
  # sink()
  flair_seg <- pred( flair_model, delta$flair,
                      gamma, alpha$flair, beta$flair,
                      lambda2$flair, a, b, m, nu2$flair, 
                      maxit$flair )
  ## Both tumor and outliers are regarded as tumor
  flair_seg$image[ flair_seg$image < -3 |
                     flair_seg$image > 0 ] <- -4L
  
  ## split t2 to grey matter and white matter
  t2_data <- splitT22( images$t2, t1ce_seg, flair_seg )
  m <- t2_data$m
  t2_model <- initEst( t2_data$label, t2_data$t2 )
  ## estimate parameters of t2 images without tumor 
  ## and CSF & necrosis
  # sink( '/home/hzhang/Documents/t2_output.txt' )
  t2_seg <- est( t2_model, delta$t2, gamma, 
                 alpha$t2[ 1 : 2 ], beta$t2[ 1 : 2 ], 
                 lambda2$t2[ 1 : 2 ],
                 m, nu2$t2, 40L )
  # sink()
  ## update beta
  sigma2 <- t2_seg$parm[ 3, ]
  beta$t2[ 1 : 2 ] <- ( alpha$t2[ 1 : 2 ] + 1 ) * sigma2
  t2_data <- split3( t2_data$t2, t2_seg, factor$t2 )
  ## update m
  m <- t2_data$m
  b <- getB( m, a )
  t2_model <- initEst( t2_data$label, t2_data$intst )
  # sink( '/home/hzhang/Documents/t2_output.txt' )
  # system.time( t2_seg <- 
  #             pred( t2_model, delta$t2, gamma, alpha$t2, beta$t2,
  #                   lambda2$t2, a, b, m, nu2$t2, maxit$t2 ) )
  # sink()
  t2_seg <- pred( t2_model, delta$t2, gamma, alpha$t2, beta$t2,
          lambda2$t2, a, b, m, nu2$t2, maxit$t2 )
  ## Get the initial results
  ## t1ce
  t1ce_image <- t1ce_seg$image
  t1ce_image[ is.na( t1ce_image ) ] <- 0L
  t1ce_image[ t1ce_image >= 1 ] <- 4L
  t1ce_image[ t1ce_image <= -4 ] <- 4L
  t1ce_image[ t1ce_image == -1 ] <- 1L
  t1ce_image[ t1ce_image == -2 ] <- 2L
  t1ce_image[ t1ce_image == -3 ] <- 3L
  ## flair
  flair_image <- flair_seg$image
  flair_image[ is.na( flair_image ) ] <- 0L
  flair_image[ flair_image >= 1 ] <- 4L
  flair_image[ flair_image <= -4 ] <- 4L
  flair_image[ flair_image == -1 ] <- 1L
  flair_image[ flair_image == -2 ] <- 2L
  flair_image[ flair_image == -3 ] <- 3L
  ## t2
  t2_image <- t2_seg$image
  t2_image[ is.na( t2_image ) ] <- 0L
  t2_image[ t2_image >= 1 ] <- 4L
  t2_image[ t2_image <= -4 ] <- 4L
  t2_image[ t2_image == -1 ] <- 1L
  t2_image[ t2_image == -2 ] <- 2L
  ## Initialize data for postprocessing
  post_data <- initPost( t1ce_image, flair_image, t2_image )
  post_seg <- postProcess( post_data )
  if( post_seg$hgg == 1 ) {
    return( post_seg$seg )
  } else {
    ## Furtherly segment edema
    further_data <- splitFthr( post_seg, t2_data )
    m <- further_data$m
    further_model <-initFther( further_data$label, further_data$intst )
    further_seg <- estF( further_model, delta$fthr, gamma,
                         alpha$fthr, beta$fthr, lambda2$fthr,
                         m, nu2$fthr, maxit$fthr )
  }
}