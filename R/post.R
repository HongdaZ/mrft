# postprocessing after segmentation
post <- function( patient, out = "SEG", infolder = "N4ITK433Z",
          ## Always four numbers for delta
          delta = 
            list( t1ce = c( -0.5, 0, 6.5, 2 ),
                  flair = c( -0.5, 0, NA_real_, 2 ),
                  t2 = c( 1, 0, NA_real_, 2 ),
                  fthr = c( 0, 0, 4, 8 ) ),
          delta_factor = 
            list( t1ce = 1,
                  flair = 1 / 1.5,
                  t2 = 1 / 3 ),
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
  ## Read segmentation results
  infile <- patient[ 1 ]
  outfile <- gsub( infolder, out, infile )
  out_t1ce_seg <- gsub( "_flair.nii.gz", "_t1ce_seg", outfile )
  t1ce_image <- readNIfTI( out_t1ce_seg, reorient = FALSE )@.Data
  out_flair_seg <- gsub( "_flair.nii.gz", "_flair_seg", outfile )
  flair_image <- readNIfTI( out_flair_seg, reorient = FALSE )@.Data
  out_t2_seg <- gsub( "_flair.nii.gz", "_t2_seg", outfile )
  t2_image <- readNIfTI( out_t2_seg, reorient = FALSE )@.Data
  ## Initialize data for postprocessing
  post_data <- initPost( t1ce_image, flair_image, t2_image )
  # sink( '/media/hzhang/ZHD-P1/result/output.txt' )
  post_seg <- postProcess( post_data, min_enh, min_enh_enc,
                           max_prop_enh_enc,
                           min_tumor, spread_add, spread_rm )
  if( sum( post_seg$image == 6, na.rm = T ) > 10 ) {
    ## With CSF inside tumor
    out_t1ce_norm <- gsub( "_flair.nii.gz", "_t1ce_norm", outfile )
    t1ce_intst <- readNIfTI( out_t1ce_norm, reorient = FALSE )@.Data
    ## Furtherly segment CSF 
    further_data <- splitFthrC( post_seg, t1ce_intst )
    m <- further_data$m
    further_model <-initFther( further_data$label, further_data$intst )
    further_seg <- est( further_model, delta$fthr[ 1 : 2 ], 
                        gamma$fthr[ 1 : 2 ],
                        alpha$fthr[ 1 : 2 ], 
                        beta$fthr[ 1 : 2  ], 
                        lambda2$fthr[ 1 : 2 ],
                        m, nu2$fthr[ 1 : 2 ], 
                        maxit$fthr )
    m <- further_seg$parm[ 2, ]
    sigma2 <- further_seg$parm[ 3, ]
    if( ( m[ 2 ] - m[ 1 ] ) /  delta$fthr[ 3 ] > 
        sqrt( sigma2[ 2 ] ) ) {
      necrosis_idx <- post_seg$image == 6 &
        further_seg$image == -2
      necrosis_idx[ is.na( necrosis_idx ) ] <- FALSE
      post_seg$image[ necrosis_idx ] <- 1L
    }
  } 
  post_seg$image[ post_seg$image == 6 ] <- NA_integer_
  post_seg$image[ post_seg$image == 5 ] <- NA_integer_
  # sink()
  if( post_seg$code != 0 ) {
    out_t2_norm <- gsub( "_flair.nii.gz", "_t2_norm", outfile )
    t2_intst <- readNIfTI( out_t2_norm, reorient = FALSE )@.Data
    ## Furtherly segment edema
    further_data <- splitFthrE( post_seg, t2_intst )
    m <- further_data$m
    further_model <-initFther( further_data$label, further_data$intst )
    further_seg <- estF( further_model, delta$fthr, gamma$fthr,
                         alpha$fthr, beta$fthr, lambda2$fthr,
                         m, nu2$fthr, maxit$fthr )
    m <- further_seg$parm[ 2, ]
    sigma2 <- further_seg$parm[ 3, ]
    if( ( m[ 2 ] - m[ 1 ] ) /  delta$fthr[ 4 ] > sqrt( sigma2[ 2 ] ) ) {
      edema_idx <- post_seg$image == 2
      edema_idx[ is.na( edema_idx ) ] <- FALSE
      post_seg$image[ edema_idx ] <- further_seg$image[ edema_idx ]
    }
  }
  # return( post_seg )
  ## Export the results to .nii images
  post_seg$image[ is.na( post_seg$image ) ] <- 0
  out_post_seg <- gsub( "_flair.nii.gz", "_post_seg", outfile )
  writeNIfTI( nifti( post_seg$image, datatype = 2 ),
              filename = out_post_seg, gzipped = TRUE )
}