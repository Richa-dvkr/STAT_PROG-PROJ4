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
  #1. If the objective or derivatives are not finite at the initial theta;
  
  try(if(is.infinite(func(theta))) stop("objective not finite at the initial theta"))
  try(if(all(is.infinite(grad(theta)))) stop("derivatives not finite at the initial theta"))
  
  ################function to get hessian matrix!!!!!!!!!!!!!!!!!!!!!!!! 
  hess=function(theta,grad)
  {
    hb1 <- grad(theta)
    Hfd <- matrix(0,2,2)       
    for (i in 1:length(theta))
    {theta_eps <- theta
    theta_eps[i] <- theta_eps[i] + eps       
    hb2 <- grad(theta_eps)                  
    Hfd[i,] <- (hb2 - hb1)/eps
    }
    return(Hfd)
  }
  ###################################################################
  
  hessian=hess(theta,grad)
  
  Hi=solve(hessian)
  
  #####################################################################
  
  
  #####function to check convergence###########################
  
  conv_grad=function(theta,grad,func,tol,fscale)
  {
    g=grad(theta)
    for(i in g)
    {
      if(abs(i)<tol*abs(func(theta)) + fscale)
      { t=TRUE}
      else
      {t=FALSE
      break}
      
    }
    return(t)
  }
  
  ###########################################################
  

  
  iter=0
  
  while(iter<=100)
  { 
    conv_test_val=conv_grad(theta,grad,func,tol,fscale)
    
    if(conv_test_val==FALSE)
    {
      obj_old=func(theta)
      gradient=grad(theta)
      Hessian=hess(theta,grad)
      
      Hi=solve(Hessian)   #######still remaining
      res <- try(chol(Hi),silent = TRUE)
      
      if(any(class(res)=="try-error"))
      {delta=-Hi%*%gradient}
      else
      {delta=Hi%*%gradient}
      
      
      new_theta=theta+delta
      print("*****")
      print(theta)
      print(delta)
      print(new_theta)
      
      obj_new=func(new_theta)
      
      print("!!!!!!!!!!!")
      print(obj_new)
      print(obj_old)
      #print(new_theta)
      #print(theta)
      
      
      no_step_half=0
      
      while(obj_old<=obj_new)
      {
        if(no_step_half<=max.half)
        {
          delta=delta/2
          new_theta=theta+delta
          obj_new=func(new_theta)
          print("*****!!!!!!")
          print(obj_new)
          print(no_step_half)
          no_step_half=no_step_half+1
          
        }
        else
        {
          print("the step fails to reduce the objective despite trying max.half step halvings")
          break
          
        }
        
        
      }
      
      if(obj_new<obj_old)
      {
        theta=new_theta
        
      }
      else{
        print("Error")
        break
        }
      
    }
    else
    {
      print(theta)
      break
    }
      
     iter=iter+1 
  }
  
  
  
}
  
  
  
  

newt(c(1,2),rb,gb)
  
  
  
  
  
  