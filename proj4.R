
#------------------------------Newton Optimizer---------------------------------

#--------------------------------Description------------------------------------

#In the following code file, the newton method's for minimization of functions
#has been implemented.For this a function named newt has been implemented. Along
#with the newt function we have also implemented functions to find the hessian
#matrix(if not passed by the user),to check whether the convergence of the function
#and a function to check whether the hessian is positive definite or not and if 
#it is not positive definite it converts it into a positive definite matrix

#--------------------------------------CODE-------------------------------------
#newt function
newt=function(theta,func,grad,hess=NULL,...,tol=1e-8,
              fscale=1,maxit=100,max.half=20,eps=1e-6)
{ 
  #Arguements
  
  # theta : vector of initial values for the optimization parameters.
  
  #func   : The objective function to minimize. Its first argument is the vector 
  #         of optimization parameters. Remaining arguments are passed from newt 
  #         using '...'
  
  #grad   : The gradient function. It has the same arguments as func but returns 
  #         the gradient vector of the objective w.r.t. the elements of parameter 
  #         vector
  
  #hess   : the Hessian matrix function. It has the same arguments as func but 
  #         returns the Hessian matrix of the objective w.r.t. the elements of 
  #         parameter vector. If not supplied then newt should obtain an 
  #         approximation to the Hessian by finite differencing of the gradient vector
  #         Default : NULL
  
  #tol    : convergence tolerance. Default = 1e-8
  
  #fscale : a rough estimate of the magnitude of func near the optimum - 
  #         used in convergence testing.Default = 1
  
  #maxit  : the maximum number of Newton iterations to try before giving up.
  #         Default = 100
  
  #maxit  : the maximum number of times a step should be halved before concluding
  #         that the step has failed to improve the objective. Default = 20
  
  #eps    : the finite difference intervals to use when a Hessian function is not provided
  #         Default: 1e-6
  
  #Returns
  # x     : a list containing
  #       f     : value of the objective function at the minimum.
  #       theta : the value of the parameters at the minimum
  #       iter  : the number of iterations taken to reach minimum
  #       g     : the gradient vector at the minimum
  #       Hi    : the inverse of the hessian matrix at the minimum
  
  
  
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
    i=1e-6     #lamda value intitiated to 1
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
      
      Hi=chol2inv(chol(Hessian))
      delta=-Hi%*%gradient  #computes delta(stepsize)
      
      new_theta=theta+delta    #computes new theta
      obj_new=func(new_theta)  #computes new obj function value
      
      no_step_half=0     #number of steps half
      while(obj_old<=obj_new || is.infinite(obj_new) )  #iterate till obj_new is 
                                                        #less than obj_old
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
      iter=iter       # number of iteration
      
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
          Hi=chol2inv(chol(Hessian))
          x=list(f=f,theta=theta,iter=iter,g=g,Hi=Hi)  # list to be return
          return(x)
        }
        
      break  #once convergence find out it breaks the loop
    } 
    
    iter=iter+1  #iter increased by 1
  }
  
}

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

newt(c(12,35),rb,gb,k=1)
newt(c(12,35),rb,gb,hb,k=1)

newt(c(1000,7),rb,gb)
newt(c(0,1),rb,gb,hb)

newt(c(20,14),rb,gb)  
newt(c(20,14),rb,gb,hb)

newt(c(1,1),rb,gb)  
newt(c(1,1),rb,gb,hb) 

newt(c(Inf,2),rb,gb)  
newt(c(1000,08),rb,gb)  

