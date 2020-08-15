# split flair to CSF & necrosis, grey matter and white matter
splitFlair3 <- function( flair, t1ce_seg ) {
  
  label <- array( -4L, dim = dim( flair ) )
  label[ ! is.nan( flair ) ] <- 0L
  
  # Find the brightest 25% percent
  q_flair <- quantile( flair, probs = .75, na.rm = T )
  bright <- flair > q_flair
  label[ bright ] <- -4L
  
  # Find CSF & necrosis
  q_flair <- quantile( flair, probs = .30, na.rm = T )
  cn <- t1ce_seg$image == -1 & flair < q_flair
  label[ cn ] <- -1L
  tbd <- label == 0
  sub_flair <- flair[ tbd ]
  sub_image <- t1ce_seg$image[ tbd ]
  
  # Find white matter
  q_flair <- quantile( sub_flair, probs =  c( .40, .50 ), na.rm = T )
  wm <- sub_image == -3L & sub_flair < q_flair[ 1 ]
  label[ tbd ][ wm ] <- -2L
  
  # Find grey matter
  gm <- sub_image == -2L & sub_flair > q_flair[ 2 ]
  label[ tbd ][ gm ] <- -3L
  
  label[ label == -4L ] <- NA_integer_
  
  # Normalize intensity
  mean_csf <- mean( flair[ label == -1 ], na.rm =  T )
  mean_white <- mean( flair[ label == -2 ], na.rm = T )
  mean_grey <- mean( flair[ label == -3 ], na.rm = T )
  
  flair <- ( flair - mean_csf ) / ( mean_grey - mean_csf )
  mean_white <- ( mean_white - mean_csf ) / ( mean_grey - mean_csf )
  res <- list( label = label, flair = flair, 
               m = c( 0, mean_white, 1 ) )
  return( res )
}
