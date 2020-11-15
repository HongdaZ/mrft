## update delta[ 3 ] for t2
updateDelta3 <- function( prop_bright,
                          t1ce_seg, flair_seg, t2_data, t2_seg ) {
  t2_brain <- t2_data$t2[ ! is.na( t2_data$t2 ) ]
  ## Find proportion of CSF & necrosis + tumor
  q_t2 <- quantile( t2_brain, probs = 1 - prop_bright, na.rm = T )
  tumor <- t2_data$t2[ t2_data$t2 > q_t2 &
                       ( t1ce_seg$image == - 4 |
                         t1ce_seg$image == - 2 ) &
                       flair_seg$image == - 4 ] 
  tumor <- tumor[ ! is.na( tumor ) ]
  x <- kmeans( tumor, 9 )$centers
  y <- order( x )
  m <- x[ y ][ 2 ]
  ## t2_seg from est
  sigma2_2 <- t2_seg$parm[ 3, 2 ]  
  m_2 <- t2_seg$parm[ 2, 2 ]
  
  delta_3 <- ( m - m_2 ) / sqrt( sigma2_2 )
  if( delta_3 > 8 ) {
    delta_3 <- 8  
  }
  return( delta_3 )
}

