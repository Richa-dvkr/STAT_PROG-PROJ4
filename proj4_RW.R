###############################################################################
#objective function
rb <- function(th,k=2) {
  k*(th[2]-th[1]**2)**2 + (1-th[1])**2
}

###############################################################################
#gradient
gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]**2),k*2*(th[2]-th[1]**2))
}

###############################################################################

newt=function(theta,func,grad,hess=NULL,...,tol=1e-8,
              fscale=1,maxit=100,max.half=20,eps=1e-6)
{
  no_iteration=1
  
  while(no_iteration<=maxit)
  {
    obj_old=func(theta)    #obj function of old point
    gradient=grad(theta)   #gradient vector of old point
    
    #A finite difference approximation to the Hessian is readily obtained by differencingthe gradient vectors
    
    hb1 <- grad(theta) ## grad at theta
    #eps from arguments
    Hfd <- matrix(0,2,2)       # finite diference Hessian
    for (i in 1:length(th0)) { ## loop over parameters
      theta_eps <- theta
      theta_eps[i] <- theta_eps[i] + eps       # increase theta_eps[i] by eps
      hb2 <- grad(theta_eps)                  #grad at th1
      Hfd[i,] <- (hb2 - hb1)/eps              #approximate second derivs
    }
    
    Hessian=Hfd   #Hessian matrix formation
    
    ###########################
    
    #need to generate function to get inverse of matrix as we cant use solve()
    
    
    
    
    Hi=solve(Hessian)
    
    
    
    
    
    ###################################
    #to check if Hessian is positive definite

    
    
    
    res <- try(chol(mat2),silent = TRUE)
    if(class(res)=="try-error")
    {
      delta=-Hi%*%gradient
    }
    else
    {
      delta=Hi%*%gradient
    }
    
    
    #####################################
    new_theta=theta+delta
    
    print(theta)
    print(delta)
    print(new_theta)
    
    obj_new=func(new_theta)
    print(obj_new)
    print(obj_old)
    

    break

    
  }
  
}



newt(c(10,.1),rb,gb)
