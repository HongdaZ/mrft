## update delta[ 3 ] for t2
updateDelta3 <- function( flair_seg, t2_seg ) {
  flair_tumor <- t2_data$t2[ flair_seg$image == -4 ]
  flair_tumor <- flair_tumor[ !is.na( flair_tumor ) ]
  ## t2_seg from est
  sigma2_2 <- t2_seg$parm[ 3, 2 ]  
  m_2 <- t2_seg$parm[ 2, 2 ]
  flair_tumor <- flair_tumor[ flair_tumor > 
                                ( m_2 + sqrt( sigma2_2 ) * 2 ) ]
  clst <- kmeans( x = flair_tumor, 
                  centers = c( min( flair_tumor ),
                               max( flair_tumor ) ) )$cluster
  m <- median( flair_tumor[ clst == 1 ] )
  
  delta_3 <- ( m - m_2 ) / sqrt( sigma2_2 ) 
  delta_3
}

