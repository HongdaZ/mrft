# Train the model

trainModel <- function( patient ) {
  lbl_intst <- splitCWG( patient, "flair" )
  
  label <- lbl_intst$label
  img <- lbl_intst$intensity
  model <- initTrn( img, label, "flair")
  
}