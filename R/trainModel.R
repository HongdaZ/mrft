# Train the model

trainModel <- function( patient, delta = 8, gamma = 1, 
                        alpha = rep( 1, 4 ), beta = rep( 1e-6, 4 ), 
                        lambda = rep( .5 / 6, 4 ), 
                        a = 5,
                        nu = rep( .5, 3 ) ) {
  
  l_intst <- splitCWG( patient )
  flair_model <- initTrn( l_intst$intensity$flair, 
                    l_intst$label$flair, "flair" )
  m <- priorMode( flair_model$info$intst, flair_model$seg )
}