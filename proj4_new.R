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
    Hfd <- matrix(0,length(theta),length(theta))           #finite diference Hessian
    for (i in 1:length(theta))     #loop over parameters
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
  
  #function check if hessian is positive finite and convert to positive
  #definite if not
  pos_func=function(mat2)
  { 
    identity_mat=diag(NCOL(mat2))  #gets identity matrix
    i=1      #lamda value intitiated to 1
    # res will gives error if matrix is not positive definite
    res <- try(chol(mat2),silent = TRUE) 
    if(any(class(res)!="try-error")) #checks if res gives error
    {
      return(mat2)   #if not it returns original matrix
    }
    else   #if it gives error--> not positive definite
    {
      while (any(class(res)=="try-error"))  #loops untill its positive definite
      {
        mat_2=mat2+i*identity_mat    #adds a multiple of the identity matrix
        i=i+1   #increase i(lambada) by 1
        res <- try(chol(mat_2),silent = TRUE)  #computes result
      }
      return(mat_2)  #returns positive definite matrix
    }
  }
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  iter=0   #intialize inter=1
  
  while(iter<=maxit)       #while loop till maxit
  { 
    #calls convergence function
    conv_test_val=conv_grad(theta,grad,func,tol,fscale,...)
    
    #*********************************************
    #checks convergence
    if(conv_test_val==FALSE)    #if convergence not reached
    { 
      if(iter==maxit)
      {stop("maxit is reached without convergence")}
      
      else
    {
      obj_old=func(theta,...)       #obj function value for theta
      gradient=grad(theta,...)      #gradient for theta
      #computes hessian using hess().using pos_func(),it checks if its positive  
      #definite or not if not it converts and returns positive definite hessian
      if(is.null(hess)){Hessian=pos_func(hess_created(theta,grad,...))}
      else{Hessian=pos_func(hess(theta,...))}
      
      
      Hi=inverse(Hessian)
      delta=-Hi%*%gradient  #computes delta(stepsize)

      new_theta=theta+delta    #computes new theta
      obj_new=func(new_theta)  #computes new obj function value
      
      no_step_half=0     #number of steps half
      
      while(obj_old<=obj_new)  #iterate till obj_new is less than obj_old
      {
        if(no_step_half<=max.half)   #checks if no_step_half exceeds max.half 
        { 
          delta=delta/2      #reducing delta by half
          new_theta=theta+delta    #new theta after reducing
          obj_new=func(new_theta,...)  #obj_new after reducing
          no_step_half=no_step_half+1   #no_step_half increased by 1
        }
        else #if no_step_half exceeds max.half we get warning
        {  
          stop("the step fails to reduce the objective despite trying
                  max.half step halvings")
        }
      }
      theta=new_theta
      
      }}
    #**************************************************
    else #if convergence reached
    { 
      g=grad(theta,...)   #optimal gradient
      f=func(theta,...)   #optimal obj function
      iter=iter       # number ofiteration
      
      #computes hessian
      if(is.null(hess)){Hessian=hess_created(theta,grad,...)}
      else{Hessian=hess(theta,...)}
      
      #computes if hessian is positive definite
      res <- try(chol(Hessian),silent = TRUE)
      #checks if res gives error(not positive definite)
      
      if(any(class(res)=="try-error"))  
        {
          stop("Hessian matrix is not finte at convergence")
          #x=list(f=f,theta=theta,iter=iter,g=g)  # list to be return
        }
      else
        {
          Hi=inverse(Hessian)
          x=list(f=f,theta=theta,iter=iter,g=g,Hi=Hi)  # list to be return
          return(x)
        }
        
      break  #once convergence find out it breaks the loop
    } 
    
    iter=iter+1  #iter increased by 1
  }
  
}



newt(c(0,2),rb,gb,k=1)

newt(c(0,2),rb,gb,hb,k=1)

newt(c(10,2),rb,gb,k=2)
newt(c(10,2),rb,gb,hb,k=2)


newt(c(28,4),rb,gb)  
newt(c(28,4),rb,gb,hb)

newt(c(1,1),rb,gb)  
newt(c(1,1),rb,gb,hb) 

newt(c(Inf,2),rb,gb)  
newt(c(1000,08),rb,gb)  




