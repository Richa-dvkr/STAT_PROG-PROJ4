
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

#########################approximation######################################


#comparing the coded gradients with finite difference approximations
fd<-th0 <- c(10,.1)
gb1 <- rb(th0)
print(gb1)
eps=1e-6
for (i in 1:length(th0)) { ## loop over parameters
  th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps
  gb2 <- rb(th1) 
  fd[i] <- (gb2 - gb1)/eps ## approximate -dl/dth[i]
}

#############

#A finite difference approximation to the Hessian is readily obtained by differencingthe gradient vectors

hb1 <- gb(th0) ## grad at th0
eps=1e-6 ## finite difference interval
Hfd <- matrix(0,2,2) ## finite diference Hessian
for (i in 1:length(th0)) { ## loop over parameters
  th1 <- th0
  th1[i] <- th1[i] + eps ## increase th0[i] by eps
  hb2 <- gb(th1) #grad at th1
  Hfd[i,] <- (hb2 - hb1)/eps ## approximate second derivs
}

Hfd



