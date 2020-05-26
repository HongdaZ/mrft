# Split CSF & necrosis, white matter and grey matter in training data
# csf = 1, wm = 2, g = 3, TBD = 0, o.w. = NA 
splitCWG <- function( patient ) {
  
  lbl <- readNifti( patient[ 2 ], internal = T )
  lbl <- as.array( lbl )
  flair <- readNifti( patient[ 1 ] )
  t2 <- readNifti( patient[ 5 ] )
  t1ce <- readNifti( patient[ 4 ] )
  
  img <- list( flair = flair, t1ce = t1ce, t2 = t2 )
  
  flair[ lbl != 0 ] <- NA
  t2[ lbl != 0 ] <- NA
  t1ce[ lbl != 0 ] <- NA
  
  
  
  # Find csf & necrosis
  q_flair <- quantile( flair, probs = c( .20, .40, .50, .01 ), na.rm = T, 
                       names = F )
  q_t2 <- quantile( t2, probs = c( .80, .01 ), na.rm = T, names = F )
  
  csf <- flair < q_flair[ 1 ] & t2 > q_t2[ 1 ]
  lbl[ csf ] <- -1L
  
  # Find white matter
  q_t1ce <- quantile( t1ce, probs = c( .60, .50, .01 ), na.rm = T, names = F )
  
  wm <- flair < q_flair[ 2 ] & t1ce > q_t1ce[ 1 ]
  lbl[ wm ] <- -2L
  
  # Find grey matter
  gm <- flair > q_flair[ 3 ] & t1ce < q_t1ce[ 2 ]
  lbl[ gm ] <- -3L
  # Remove the darkest 1% of the voxels
  lbl_flair <- lbl
  lbl_flair[ flair < q_flair[ 4 ] ] <- NA_integer_
  lbl_t1ce <- lbl
  lbl_t1ce[ t1ce < q_t1ce[ 3 ] ] <- NA_integer_
  lbl_t2 <- lbl
  lbl_t2[ t2 < q_t2[ 2 ] ] <- NA_integer_
  
  # storage.mode( lbl_flair ) <- "integer"
  # storage.mode( lbl_t1ce ) <- "integer"
  # storage.mode( lbl_t2 ) <- "integer"
  
  lbl <- list( flair = lbl_flair, t1ce = lbl_t1ce, t2 = lbl_t2 )
  list( label = lbl, intensity = img )
}