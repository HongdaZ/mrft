# estimate parameters of t1ce or flair images without tumor
est3 <- function( model, delta, gamma, 
                  alpha, beta, lambda, m, nu2, maxit ) {
  .Call( "est3", model, delta, gamma, 
         alpha, beta, lambda, m, nu2, maxit )
}