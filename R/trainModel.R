# Train the model

trainModel <- function( patient, delta = 8, gamma = 1, 
                        alpha = rep( 1, 4 ), beta = rep( 1e-6, 4 ), 
                        lambda2 = rep( .25 / 72, 4 ), 
                        a = 5,
                        nu2 = rep( .25, 3 ) ) {
  
  l_intst <- splitCWG( patient ) # splitCWGX for testing
  flair_model <- initTrn( l_intst$intensity$flair, 
                    l_intst$label$flair, "flair" )
  # -1, -2, -3, > 0
  m <- priorMode( flair_model$info$intst, flair_model$seg )
  b <- getB( m, a )
  par_flair <- estParm( flair_model, delta, gamma, 
                        alpha, beta, lambda2, a, b, m, nu2 )
}
