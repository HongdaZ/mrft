# split T1ce to CSF & necrosis, grey matter and white matter
splitT1ce3 <- function( t1ce, flair ) {
  
  label <- array( -4L, dim = dim( flair ) )
  label[ ! is.nan( t1ce ) ] <- 0L
  q_flair <- quantile( flair, probs = .15, na.rm = T )
  q_t1ce <- quantile( t1ce, probs = .15, na.rm = T )
  
  # find CSF & necrosis
  label[ flair < q_flair & t1ce < q_t1ce ] <- -1L
  # Remove bright regions
  q_t1ce <- quantile( t1ce, probs = .90, na.rm = T )
  bright <- t1ce > q_t1ce
  label[ bright ] <- -4L
  
  tbd <- label == 0
  
  sub_flair <- flair[ tbd ]
  sub_t1ce <- t1ce[ tbd ]
  
  q_flair <- quantile( sub_flair, probs = c( .40, .50 ), 
                       na.rm =  T )
  q_t1ce <- quantile( sub_t1ce, probs =  c( .50, .60 ),
                      na.rm =  T )
  # find white matter
  label[ tbd ][ sub_flair < q_flair[ 1 ] & sub_t1ce > q_t1ce[ 2 ] ] <- -3L
  # find grey matter
  label[ tbd ][ sub_flair > q_flair[ 2 ] & sub_t1ce < q_t1ce[ 1 ] ] <- -2L
  
  label[ label == -4L ] <- NA_integer_
  
  # normalize intensity
  mean_csf <- mean( t1ce[ label == -1 ], na.rm = T )
  mean_grey <- mean( t1ce[ label == -2 ], na.rm = T )
  mean_white <- mean( t1ce[ label == - 3 ], na.rm = T )
  
  t1ce <- ( t1ce - mean_csf ) / ( mean_white - mean_csf ) 
  mean_grey <- ( mean_grey - mean_csf ) / ( mean_white - mean_csf ) 
  res <- list( label = label, t1ce = t1ce,
               m = c( 0, mean_grey, 1 ) )
  return( res )
}
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
# split t1ce or flair images for prediction
split4 <- function( x, x_seg, x_factor ) {
  label <- array( -4L, dim = dim( x ) )
  label[ ! is.nan( x ) ] <- 0L
  label[ x_seg$image == -1 ] <- -1L
  label[ x_seg$image == -2 ] <- -2L
  label[ x_seg$image == -3 ] <- -3L
  m_3 <- x_seg$parm[ 2, 1 ]
  sigma2_3 <- x_seg$parm[ 3, 1 ]
  # threshold for enhancing tumor core
  m_4 <- m_3 + sqrt( sigma2_3 ) * x_factor
  label[ x > m_4 ] <- 4L
  m_1 <- x_seg$parm[ 2, 3 ]
  m_2 <- x_seg$parm[ 2, 2 ]
  
  label[ label == -4L ] <- NA_integer_
  
  res <- list( label = label, intst = x, 
               m = c( m_1, m_2, m_3, m_4 ) )
  return( res )
}