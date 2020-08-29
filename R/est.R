## estimate parameters of t1ce or flair images without tumor
est3 <- function( model, delta, gamma, 
                  alpha, beta, lambda2, m, nu2, maxit ) {
  .Call( "est3", model, delta, gamma, 
         alpha, beta, lambda2, m, nu2, maxit )
}

## estimate parameters of t2 images without tumor and CSF & necrosis
est2 <- function( model, delta, gamma, 
                  alpha, beta, lambda2, m, nu2, maxit ) {
  .Call( "est2", model, delta, gamma, 
         alpha, beta, lambda2, m, nu2, maxit )
}