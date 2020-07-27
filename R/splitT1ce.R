# split T1ce to CSF & necrosis, grey matter and white matter
splitT1ce <- function( t1ce, flair ) {
  
  label <- array( -4L, dim = dim( t1ce ) )
  label[ ! is.nan( t1ce ) ] <- 0L
  q_flair <- quantile( flair, probs = .30, na.rm = T )
  q_t1ce <- quantile( t1ce, probs = .30, na.rm = T )
  
  # find CSF & necrosis
  label[ flair < q_flair & t1ce < q_t1ce ] <- -1L
  # Remove bright regions
  q_t1ce <- quantile( t1ce, probs = .90, na.rm = T )
  bright <- t1ce > q_t1ce
  t1ce[ bright ] <- NaN
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
  t1ce[ is.na( label ) ] <- NaN
  res <- list( label = label, t1ce = t1ce )
  return( res )
}