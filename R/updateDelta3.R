## update delta[ 3 ] for t2
updateDelta3 <- function( flair_seg, t2_data, t2_seg ) {
  flair_tumor <- t2_data$t2[ flair_seg$image == - 4 ]
  flair_tumor <- flair_tumor[ !is.na( flair_tumor ) ]
  ## t2_seg from est
  sigma2_2 <- t2_seg$parm[ 3, 2 ]  
  m_2 <- t2_seg$parm[ 2, 2 ]
  flair_tumor <- flair_tumor[ flair_tumor > 
                                ( m_2 + sqrt( sigma2_2 ) * 4 ) ]
  clst <- kmeans( x = flair_tumor, 
                  centers = summary( flair_tumor )[ - 4 ] )$cluster

  m <- tapply( flair_tumor, clst, median )
  
  delta_3 <- ( m - m_2 ) / sqrt( sigma2_2 )
  delta_3 <- kmeans( x = delta_3,
                   c( min( delta_3 ),
                      median( delta_3 ),
                      max( delta_3 ) ) )
  delta_3 <- delta_3$centers[ 2 ]
  delta_3
}

