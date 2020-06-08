# fit model to estimate mu1, mu2, mu3 and sigma's
estParm <- function( model, delta, gamma, 
                            alpha, beta, lambda, a, b, m, nu ) {
  .Call( "estParm", model, delta, gamma, 
         alpha, beta, lambda, a, b, m, nu )
}