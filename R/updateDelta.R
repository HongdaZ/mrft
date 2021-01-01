## update delta[ 3 ] for t2
updateDelta3T2 <- function( prop_bright,
                          t1ce_image, flair_image, t2_data, t2_seg ) {
  t2_brain <- t2_data$t2[ ! is.na( t2_data$t2 ) ]
  ## Find proportion of CSF & necrosis + tumor
  q_t2 <- quantile( t2_brain, probs = 1 - prop_bright, na.rm = T )
  tumor <- t2_data$t2[ t2_data$t2 > q_t2 &
                       ( t1ce_image == 4 |
                         t1ce_image == 2 ) &
                       flair_image == 4 ] 
  tumor <- tumor[ ! is.na( tumor ) ]
  start <- quantile( tumor, probs = seq( 0, 1, length.out = 9 ) )
  x <- kmeans( tumor, start )$centers
  y <- order( x )
  m <- x[ y ][ 2 ]
  ## t2_seg from est
  sigma2_2 <- t2_seg$parm[ 3, 2 ]  
  m_2 <- t2_seg$parm[ 2, 2 ]
  
  delta_3 <- ( m - m_2 ) / sqrt( sigma2_2 )
  if( delta_3 > 12.5 ) {
    delta_3 <- 12.5  
  }
  return( delta_3 )
}

## update delta[ 3 ] for flair
updateDelta3Flair <- function( t1ce_image, flair_data, flair_seg ) {
  flair_brain <- flair_data$flair[ ! is.na( flair_data$flair ) ]
  ## Find grey matter + tumor
  q_flair <- quantile( flair_brain, probs = 0.70, na.rm = T )
  tumor <- flair_data$flair[ flair_data$flair > q_flair &
                             ( t1ce_image == 4 |
                               t1ce_image == 2 ) ] 
  tumor <- tumor[ ! is.na( tumor ) ]
  start <- quantile( tumor, probs = seq( 0, 1 - 0.12, length.out = 8 ) )
  x <- kmeans( tumor, start )$centers
  y <- order( x )
  m <- x[ y ][ 3 ]
  ## t2_seg from est
  sigma2_3 <- flair_seg$parm[ 3, 3 ]  
  m_3 <- flair_seg$parm[ 2, 3 ]
  
  delta_3 <- ( m - m_3 ) / sqrt( sigma2_3 )
  if( delta_3 > 13 ) {
    delta_3 <- 13  
  } else if( delta_3 < 4 ) {
    delta_3 <- 4
  }
  return( delta_3 )
}

## Update delta[ 1, 2 ] or delta[ 1 ] for Flair or T2
updateDelta12 <- function( m_u, sigma_u, m_l, sigma_l, shift_t ) {
  res <- ( ( shift_t / 2 ) ^ 2 - 
             ( ( m_u + shift_t / 2 * sqrt( sigma_u ) - m_l ) /
                 sqrt( sigma_l ) ) ^ 2 ) / 2  
  if( res < 0 ) {
    res <- 0
  }
  res
}
