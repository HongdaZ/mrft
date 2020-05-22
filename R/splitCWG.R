# Split CSF & necrosis, white matter and grey matter in training data
# csf = 1, wm = 2, g = 3, TBD = 0, o.w. = NA 
splitCWG <- function( patient ) {
  
  lbl <- readNifti( patient[ 2 ] )
  
  flair <- readNifti( patient[ 1 ] )
  flair[ lbl != 0 ] <- NA
  
  t2 <- readNifti( patient[ 5 ] )
  t2[ lbl != 0 ] <- NA
  
  t1ce <- readNifti( patient[ 4 ] )
  t1ce[ lbl != 0 ] <- NA
  
  # Find csf & necrosis
  q_flair <- quantile( flair, probs = .20, na.rm = T )
  q_t2 <- quantile( t2, probs = .80, na.rm = T )
  
  csf <- flair < q_flair & t2 > q_t2
  lbl[ csf ] <- -1
  
  # Find white matter
  q_flair <- quantile( flair, probs = .40, na.rm = T )
  q_t1ce <- quantile( t1ce, probs = .60, na.rm = T )
  
  wm <- flair < q_flair & t1ce > q_t1ce
  lbl[ wm ] <- -2
  
  # Find grey matter
  q_flair <- quantile( flair, probs = .50, na.rm = T )
  q_t1ce <- quantile( t1ce, probs = .50, na.rm = T )
  
  gm <- flair > q_flair & t1ce < q_t1ce
  lbl[ gm ] <- -3
  lbl
}