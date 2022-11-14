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
  
  #function check and convert hessian to positive definite
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
    #computes hessian using hess().using pos_func(),it checks if its positive  
    #definite or not if not it converts and returns positive definite hessian
      Hessian=pos_func(hess(theta,grad))
      
      Hi=chol2inv(chol(Hessian))
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
      g=grad(theta)   #optimal gradient
      f=func(theta)   #optimal obj function
      iter=iter       # number ofiteration
      
      Hessian=hess(theta,grad) #computes hessian
      #computes if hessian is positive definite
      res <- try(chol(Hessian),silent = TRUE)
      #checks if res gives error(not positive definite)
      if(any(class(res)=="try-error"))  
      {
        warning("Hessian matrix is not finte at convergence")
        x=list(f=f,theta=theta,iter=iter,g=g)  # list to be return
      }
      else
      {
        Hi=chol2inv(chol(Hessian))
        x=list(f=f,theta=theta,iter=iter,g=g,Hi=Hi)  # list to be return
      }
      
      return(x)
      
      break  #once convergence find out it breaks the loop
    } 
    
    iter=iter+1  #iter increased by 1
  }
  
  
  
  }
      
      
      
newt(c(10,.1),rb,gb)
newt(c(1,1),rb,gb)
  
  
  
  
  
  