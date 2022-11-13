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
#newt function

newt=function(theta,func,grad,hess=NULL,...,tol=1e-8,
              fscale=1,maxit=100,max.half=20,eps=1e-6)
{
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #checks if the objective are not finite at the initial theta;
  try(if(is.infinite(func(theta)))stop("objective not finite 
                                       at the initial theta"))
  
  #checks if the derivatives are not finite at the initial theta;
  try(if(all(is.infinite(grad(theta))))stop("derivatives not 
                                            finite at the initial theta"))
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  #function to form Hessian Matrix
  hess=function(theta,grad)
  {
    hb1 <- grad(theta)             #grad at theta
    Hfd <- matrix(0,2,2)           #finite diference Hessian
    for (i in 1:length(theta))     #loop over parameters
    {theta_eps <- theta
    theta_eps[i] <- theta_eps[i] + eps   #increase th0[i] by eps      
    hb2 <- grad(theta_eps)               #grad at theta
    Hfd[i,] <- (hb2 - hb1)/eps           #approximate second derivatives
    }
    return(Hfd)      #returns hessian matrix
  }
  
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  #function to form inverse of matrix
  
  #function to find minors
  minor <- function(mat){
    rows <- dim(mat)[1]  #number of rows of matrix  
    cols <- dim(mat)[2]  #number of columns of matrix 
    minor <- matrix(0,nrow=rows,ncol=cols)      #creating empty minor matrix
    for (r in range(1,rows)){      #looping over rows
      for (c in range(1,cols)){    #looping over columns
        minor[r,c] = mat[-r,-c]}}  #finding minors
    return(minor)    #returns minor
  }
  
  #function to find cofactors
  cofactor <- function(mat){
    rows <- dim(mat)[1]  #number of rows of matrix 
    cols <- dim(mat)[2]  #number of columns of matrix
    cofac <- matrix(0,nrow=rows,ncol=cols)      #creating empty cofactors matrix
    for (r in range(1,rows)){         #looping over rows
      for (c in range(1,cols)){       #looping over columns
        e <- (r+c)          #sum of index
        cofac[r,c] = ((-1) ** e) * minor(mat)[r,c]}}  #finds cofactors
    return(cofac)   #returns cofactors
  }
  
  #function to find inverse
  inverse <- function(mat){
    cofac_t <- t(cofactor(mat))   #cofactor transpose
    inv <- cofac_t /det(mat)      #inverse matrix formula
    return(inv)     #returns inverse
  }
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  #function to check convergence
  conv_grad=function(theta,grad,func,tol,fscale)
  {
    gradient=grad(theta)   #forms gradient
    for(ele in gradient)   #loop over elements of gradient
    {
      if(abs(ele)<tol*abs(func(theta)) + fscale)  #checks convergence
      { result=TRUE}     
      else
      {result=FALSE
      break}
      
    }
    return(result)    #returns result
  }
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  iter=1   #intialize inter=1
  
  while(iter<=maxit)       #while loop till maxit
  { 
    #calls convergence function
    conv_test_val=conv_grad(theta,grad,func,tol,fscale)
    
    #********************************************
    #checks if maxit reach without convergence
    if(iter==maxit && conv_test_val==FALSE)
    {
      warning("maxit is reached without convergence")
      break
    }
    
    #*********************************************
    #checks convergence
    if(conv_test_val==FALSE)    #if convergence not reached
    {
      obj_old=func(theta)       #obj function value for theta
      gradient=grad(theta)      #gradient for theta
      Hessian=hess(theta,grad)  #hessian for theta
      Hi=inverse(Hessian)       #inverse of hessian matrix
      
      res = try(chol(Hi),silent = TRUE) #computes if hessian is positive definite
      
      
      #checks if res gives error(not positive definite)
      if(any(class(res)=="try-error")) 
      {delta=-Hi%*%gradient}    #computes delta(stepsize) 
      else
      {delta=Hi%*%gradient}    #computes delta(stepsize)
      
      
      new_theta=theta+delta    #computes new theta
      obj_new=func(new_theta)  #computes new obj function value
      
      
      no_step_half=0     #number of steps half
      
      while(obj_old<=obj_new)  #iterate till obj_new is less than obj_old
      {
        if(no_step_half<=max.half)   #checks if no_step_half exceeds max.half 
        {
          delta=delta/2      #reducing delta by half
          new_theta=theta+delta    #new theta after reducing
          obj_new=func(new_theta)  #obj_new after reducing
          no_step_half=no_step_half+1   #no_step_half increased by 1
        }
        else #if no_step_half exceeds max.half we get warning
        {  
          warning("the step fails to reduce the objective despite trying
                  max.half step halvings")
          break   #terminates loop
        }
  
      }
      
      if(obj_new<obj_old)  #if obj_new less than obj_old without any errors
      { 
        theta=new_theta
        }
      else{
        break
        }
      
    }
    #**************************************************
    else #if convergence reached
    {
      Hessian=hess(theta,grad) #computes hessian
      Hi=inverse(Hessian)       #inverse of hessian matrix
      #computes if hessian is positive definite
      res <- try(chol(Hi),silent = TRUE)
      #checks if res gives error(not positive definite)
      if(any(class(res)=="try-error"))  
      {
        warning("Hessian matrix is not finte at convergence")
        Hi=NA   #initialize Hi with NA
      }
      
      #!!!!!!!!!!!!!!!!!!!!return list!!!!!!!!!!!!!!!!!
      g=grad(theta)   #optimal gradient
      f=func(theta)   #optimal obj function
      iter=iter       # number ofiteration
      x=list(f=f,theta=theta,iter=iter,g=g,Hi=Hi)  # list to be return
      return(x)

      break  #once convergence find out it breaks the loop
    } 

    iter=iter+1  #iter increased by 1
  }
  
  
  
}

newt(c(10,.1),rb,gb)
newt(c(1,5),rb,gb)

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  



















