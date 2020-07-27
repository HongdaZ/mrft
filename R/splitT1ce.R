# split T1ce to CSF & necrosis, grey matter and white matter
splitT1ce <- function( t1ce, flair ) {
  
  label <- array( -4L, dim = dim( t1ce ) )
  label[ ! is.nan( t1ce ) ] <- 0L
  q_flair <- quantile( flair, probs = .30, na.rm = T )
  q_t1ce <- quantile( t1ce, probs = .30, na.rm = T )
  
  # find CSF & necrosis
  label[ flair < q_flair & t1ce < q_t1ce ] <- -1L
  # Remove bright regions
  q_t1ce <- quantile( t1ce, probs = .80, na.rm = T )
  label[ t1ce > q_t1ce ] <- -4L
  
  tbd <- label == 0
  q_flair <- quantile( flair[ tbd ], probs = c( .40, .50 ), 
                       na.rm =  T )
  q_t1ce <- quantile( t1ce[ tbd ], probs =  c( .50, .60 ),
                      na.rm =  T )
  # find white matter
  label[ flair[ tbd ] < q_flair[ 1 ] & t1ce[ tbd ] > q_t1ce[ 2 ] ] <- -3L
  # find grey matter
  label[ flair[ tbd ] > q_flair[ 2 ] & t1ce[ tbd ] < q_t1ce[ 1 ] ] <- -2L
  
  label[ label == -4L ] <- NA_integer_
  t1ce[ is.na( label ) ] <- NaN
  res <- list( label = label, t1ce = t1ce )
  return( res )
}