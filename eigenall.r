eigenall<-function(A, zero=TRUE) # A is a matrix
{
    ev <- eigen(A)
    # R sorts eigenvalues in decreasing order, according to Mod(values)
    #  ususally dominant eigenvalue is first (ev$values[1]), except for imprimitive matrices with d eigenvalues of equal modulus
    # this should work for most cases
    lmax <- which.max(Re(ev$values))
    lambda <- Re(ev$values[lmax])
    ## Damping ratio. Use second eigenvalue
    # dr<- lambda/abs(ev$values[2])
     ## OR  second largest magnitude in case of ties using rle - round needed for imprimitive matrices
    dr<-rle(round(Mod(ev$values), 5  ))$values
    dr<-dr[1]/dr[2]
    
    W <- ev$vectors
    w <- abs(Re(W[, lmax]))
    ## check if matrix is singular-and output NAs rather than stop (better for loops and bootstrapping)
    V <- try(Conj(solve(W)), silent=TRUE)
    if (class(V) == "try-error") {
      eigenall <- list(lambda1 = lambda, stable.stage = w/sum(w), 
        sensitivities = A*NA, elasticities = A*NA, repro.value = w*NA, 
        damping.ratio = dr)
                           }
    else{ 
    v <- abs(Re(V[lmax, ]))
    S=(cbind(v/v[1]) %*% w/sum(w))/(((v/v[1]) %*%w/sum(w))[1]) # sensitivity matrix
    E=A*S/lambda #The elasticity matrix


    eigenall <- list(lambda1 = lambda, stable.stage = w/sum(w), 
         repro.value = v/v[1], 
        damping.ratio = dr,sensitivities =S,elasticities =E )
  }
    
}
