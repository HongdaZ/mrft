# predict the labels for T1ce or T2 images
pred4 <- function( model, delta, gamma,
                  alpha, beta, lambda2, a, b, m, nu2, maxit ) {
  .Call( "pred4", model, delta, gamma, 
         alpha, beta, lambda2, a, b, m, nu2, maxit )
}