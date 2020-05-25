# Train the model

trainModel <- function( patient ) {
  
  l_intst <- splitCWG( patient )
  flair_model <- initTrn( l_intst$intensity$flair, 
                    l_intst$label$flair, "flair")
  
}