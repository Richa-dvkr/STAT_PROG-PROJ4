mat2.data <- c(-4,1,2,1)
mat2 <- matrix(mat2.data,nrow=2)
mat2


pos_func=function(mat2)
{ 
  identity_mat=diag(NCOL(mat2))  #gets identity matrix
  i=1      #lamda value intitiated to 1
  res <- try(chol(mat2),silent = TRUE) 
  if(any(class(res)!="try-error"))
  {
    return(mat2)
  }
  else
  {
    while (any(class(res)=="try-error"))
    {
      mat_2=mat2+i*identity_mat
      i=i+1
      res <- try(chol(mat_2),silent = TRUE)
    }
    return(mat_2)
    
  }
  
}
pos_func(mat2)

########################

#objective function
rb <- function(th,k=2) {
  k*(th[2]-th[1]**2)**2 + (1-th[1])**2
}


#gradient
gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]**2),k*2*(th[2]-th[1]**2))
}

#hessian
hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]**2) - 4*th[1]**2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}
#######################





hh

hh=hb(c(10,.1))
solve(hh)

pos_func=function(mat2)
{ t=diag(NCOL(mat2))
  i=1
  res <- try(chol(mat2),silent = TRUE)
  
  while (any(class(res)=="try-error")) {
    mat_2=mat2+i*t
    i=i+1
    res <- try(chol(mat_2),silent = TRUE)
    
    }
  return(mat_2)
}
pos_func(hh)

