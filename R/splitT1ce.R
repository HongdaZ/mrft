# split T1ce to CSF & necrosis, grey matter and white matter
splitT1ce3 <- function( t1ce, flair ) {
  
  label <- array( -4L, dim = dim( t1ce ) )
  label[ ! is.nan( t1ce ) ] <- 0L
  q_flair <- quantile( flair, probs = .30, na.rm = T )
  q_t1ce <- quantile( t1ce, probs = .30, na.rm = T )
  
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
# split t1ce images for prediction
splitT1ce4 <- function( t1ce, t1ce_seg ) {
  label <- array( NA_integer_, dim = dim( t1ce ) )
  label[ ! is.nan( t1ce ) ] <- 0L
  label[ t1ce_seg$image == -1 ] <- -1L
  label[ t1ce_seg$image == -2 ] <- -2L
  label[ t1ce_seg$image == -3 ] <- -3L
  m_3 <- t1ce_seg$parm[ 2, 1 ]
  sigma2_3 <- t1ce_seg$parm[ 3, 1 ]
  # threshold for enhancing tumor core
  m_4 <- m_3 + sqrt( sigma2_3 ) * 12
  label[ t1ce > m_4 ] <- 4L
  m_1 <- t1ce_seg$parm[ 2, 3 ]
  m_2 <- t1ce_seg$parm[ 2, 2 ]
  
  res <- list( label = label, t1ce = t1ce, 
               m = c( m_1, m_2, m_3, m_4 ) )
  return( res )
}