# segment the MR images 
# beta and nu2 influenced by normalization

segment <- function( patient, out = "SEG", infolder = "N4ITK433Z",
                     ## Always four numbers for delta
                     delta = 
                       list( t1ce = c( -3, -2, 6, 6 ),
                             flair = c( 0, 0, 8, 8 ),
                             t2 = c( 3, 0, NA_real_, NA_real_ ),
                             fthr = c( 0, 0, 8, 0 ) ),
                     gamma = 0.80,
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
  res <- matrix( nrow = 2, ncol = 3 )
  row.names( res ) <- c( "m", "sigma2" )
  colnames( res ) <- c( "t1ce", "flair", "t2" )
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
  res[ , 1 ] <- t1ce_seg$parm[ c( 2, 3 ), 3 ]
  ## update beta
  sigma2 <- t1ce_seg$parm[ 3, ]
  beta$t1ce[ 1 : 3 ] <- ( alpha$t1ce[ 1 : 3 ] + 1 ) * ( sigma2 )
  t1ce_data <- split4( t1ce_data$t1ce, t1ce_seg, delta$t1ce[ 3 ] )
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
  res[ , 2 ] <- flair_seg$parm[ c( 2, 3 ), 3 ]
  ## update beta
  sigma2 <- flair_seg$parm[ 3, ]
  beta$flair[ 1 : 3 ] <- ( alpha$flair[ 1 : 3 ] + 1 ) * ( sigma2 )
  flair_data <- split4( flair_data$flair, flair_seg, delta$flair[ 3 ] )
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
  prop_bright <- propBright( t1ce_seg, flair_seg, images$t2 )
  t2_data <- splitT22( prop_bright, images$t2, t1ce_seg, flair_seg )
  m <- t2_data$m
  t2_model <- initEst( t2_data$label, t2_data$t2 )
  ## estimate parameters of t2 images without tumor 
  ## and CSF & necrosis
  # sink( '/home/hzhang/Documents/t2_output.txt' )
  t2_seg <- est( t2_model, delta$t2, gamma, 
                 alpha$t2[ 1 : 2 ], beta$t2[ 1 : 2 ], 
                 lambda2$t2[ 1 : 2 ],
                 m, nu2$t2, 40L )
  res[ , 3 ] <- t2_seg$parm[ c( 2, 3 ), 2 ]
  # sink()
  ## update delta
  if( is.na( delta$t2[ 3 ] ) ) {
    delta$t2[ c( 3, 4 ) ] <- updateDelta3( prop_bright,
                                           t1ce_seg, flair_seg,
                                           t2_data, t2_seg ) 
  }
  ## update beta
  sigma2 <- t2_seg$parm[ 3, ]
  beta$t2[ 1 : 2 ] <- ( alpha$t2[ 1 : 2 ] + 1 ) * sigma2
  t2_data <- split3( t2_data$t2, t2_seg, delta$t2[ 3 ] )
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
  # sink( '/media/hzhang/ZHD-P1/result/output.txt' )
  post_seg <- postProcess( post_data, min_enh, max_prop_enh_enc,
                           min_tumor, spread, 
                           min_prop_tumor_nbr )
  # sink()
  if( post_seg$code != 0 ) {
    ## Furtherly segment edema
    further_data <- splitFthr( post_seg, t2_data )
    m <- further_data$m
    further_model <-initFther( further_data$label, further_data$intst )
    further_seg <- estF( further_model, delta$fthr, gamma,
                         alpha$fthr, beta$fthr, lambda2$fthr,
                         m, nu2$fthr, maxit$fthr )
    m <- further_seg$parm[ 2, ]
    sigma2 <- further_seg$parm[ 3, ]
    if( ( m[ 2 ] - m[ 1 ] ) /  delta$fthr[ 3 ] > sqrt( sigma2[ 2 ] ) ) {
      edema_idx <- post_seg$image == 2
      edema_idx[ is.na( edema_idx ) ] <- FALSE
      post_seg$image[ edema_idx ] <- further_seg$image[ edema_idx ]
    }
  }
  # return( post_seg )
  ## Export the results to .nii images
  infile <- patient[ 1 ]
  outfile <- gsub( infolder, out, infile )
  out_t1ce <- gsub( "_flair.nii.gz", "_t1ce_seg", outfile )
  writeNIfTI( nifti( t1ce_image, datatype = 2 ),
              filename = out_t1ce, gzipped = TRUE )
  out_flair <- gsub( "_flair.nii.gz", "_flair_seg", outfile )
  writeNIfTI( nifti( flair_image, datatype = 2 ),
              filename = out_flair, gzipped = TRUE )
  out_t2 <- gsub( "_flair.nii.gz", "_t2_seg", outfile )
  writeNIfTI( nifti( t2_image, datatype = 2 ),
              filename = out_t2, gzipped = TRUE )
  post_seg$image[ is.na( post_seg$image ) ] <- 0
  out_post <- gsub( "_flair.nii.gz", "_post_seg", outfile )
  writeNIfTI( nifti( post_seg$image, datatype = 2 ),
              filename = out_post, gzipped = TRUE )
  ## Export normalized images
  t1ce_intst <- t1ce_data$intst
  out_t1ce <- gsub( "_flair.nii.gz", "_t1ce_norm", outfile )
  writeNIfTI( nifti( t1ce_intst, datatype = 64 ) ,
              out_t1ce,
              gzipped = T )
  flair_intst <- flair_data$intst
  out_flair <- gsub( "_flair.nii.gz", "_flair_norm", outfile )
  writeNIfTI( nifti( flair_intst, datatype = 64 ) ,
              out_flair,
              gzipped = T )
  t2_intst <- t2_data$intst
  out_t2 <- gsub( "_flair.nii.gz", "_t2_norm", outfile )
  writeNIfTI( nifti( t2_intst, datatype = 64 ) ,
              out_t2,
              gzipped = T )
  data_file <- gsub( "_flair.nii.gz", ".RData", outfile )
  save( res, file = data_file )
}