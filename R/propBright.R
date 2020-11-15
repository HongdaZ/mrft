## Find a rough proportion of bright regions in T2
propBright <- function( t1ce_seg, flair_seg ) {
  n_cnt <- sum( t1ce_seg$image == - 1 | 
                  t1ce_seg$image == - 4 |
                  flair_seg$image == - 1 |
                  flair_seg$image == - 4, na.rm = T  ) 
  n_b <- sum( ! is.na( t2_data$t2 ) )
  prop_bright <- n_cnt / n_b 
  return( prop_bright )
}