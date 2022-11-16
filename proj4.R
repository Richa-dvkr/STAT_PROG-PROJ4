###############################################################################
#objective function

rb <- function(th,k=2) {
  k*(th[2]-th[1]**2)**2 + (1-th[1])**2
}

###############################################################################
###############################################################################
#gradient function

gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]**2),k*2*(th[2]-th[1]**2))
}

###############################################################################
###############################################################################
#hessian matrix
hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]**2) - 4*th[1]**2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}


###############################################################################
###############################################################################
#newt function

newt=function(theta,func,grad,hess=NULL,...,tol=1e-8,
              fscale=1,maxit=100,max.half=20,eps=1e-6)
{ 
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #checks if the objective are not finite at the initial theta;
  if(is.infinite(func(theta,...)))
  {
    stop("objective not finite at the initial theta")
  }
  
  #checks if the derivatives are not finite at the initial theta;
  if(all(is.infinite(grad(theta,...))))
  {
    stop("derivatives not finite at the initial theta")
  }
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  #function to form Hessian Matrix
  hess_created=function(theta,grad,...)
  {
    hb1 <- grad(theta,...)             #grad at theta
    Hfd <- matrix(0,length(theta),length(theta))     #finite difference Hessian
    for (i in 1:length(theta))      #loop over parameters
    {theta_eps <- theta
    theta_eps[i] <- theta_eps[i] + eps   #increase th0[i] by eps      
    hb2 <- grad(theta_eps,...)               #grad at theta
    Hfd[i,] <- (hb2 - hb1)/eps           #approximate second derivatives
    }
    return(Hfd)      #returns hessian matrix
  }
  
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  #function to check convergence
  conv_grad=function(theta,grad,func,tol,fscale,...)
  {
    gradient=grad(theta,...)   #forms gradient
    for(ele in gradient)   #loop over elements of gradient
    {
      if(abs(ele)< ((fscale + abs(func(theta,...)))*tol))  #checks convergence
      { result=TRUE}     
      else
      {result=FALSE
      break}
      
    }
    return(result)    #returns result
  }
  
  pos_func=function(mat2)
  { 
    identity_mat = diag(NCOL(mat2))  #gets identity matrix
    i=eps     #lamda value intitiated to 1
    # res will gives error if matrix is not positive definite
    res <- try(chol(mat2),silent = TRUE) 
    if(any(class(res)!="try-error")) #checks if res gives error
    {
      return(mat2)   #if not it returns original matrix
    }
    else   #if it gives error--> not positive definite
    {
      while (any(class(res)=="try-error"))  #loops until its positive definite
      { 
        mat2=(mat2+t(mat2))/2
        mat_2=mat2+i*identity_mat    #adds a multiple of the identity matrix
        i=i*10                       #increase i(lamda) by 1
        res <- try(chol(mat_2),silent = TRUE)  #computes result
      }
      return(mat_2)  #returns positive definite matrix
    }
  }
  
}  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!