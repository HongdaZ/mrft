## Initialize data for training
initTrn <- function( label_, intst_, modality ) {

  if( modality == "flair" ) {
    # label[ label == 1 | label == 4 ] <- NA_integer_
    label <- changeD( label_, 1L, 4L )
  } else if( modality == "t1ce" ) {
    # label[ label == 1 | label == 2 ] <- NA_integer_
    label <- changeD( label_, 1L, 2L )
  } else if( modality == "t2" ) {
    # label[ label == 4 ] <- NA_integer_
    label <- changeD( label_, 4L, NA_integer_ )
  }
  info <- indexMat( intst_, label )
  seg <- label[ info$idx ]
  padding <- rep( 0L, length( seg ) )
  seg_mat <- rbind( seg, padding )
  
  res <- list( info = info, seg = seg_mat )
  return( res )
}