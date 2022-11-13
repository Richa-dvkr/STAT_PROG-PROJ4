#objective function
Dtheta <- function(th,k=2) {
  k*(th[2]-th[1]**2)**2 + (1-th[1])**2
}

#gradient 
delta_Dtheta <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]**2),k*2*(th[2]-th[1]**2))
}

#hessian matrix function
deltasq_Dtheta <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]**2) - 4*th[1]**2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}


###############################################################

minor <- function(mat){
  rows <- dim(mat)[1]
  cols <- dim(mat)[2]
  minor <- matrix( 0 , nrow = rows , ncol = cols )
  for (r in range(1,rows)){
    for (c in range(1,cols)){
      minor[r,c] = mat[-r,-c]}}
  return(minor)
}


cofactor <- function(mat){
  rows <- dim(mat)[1]
  cols <- dim(mat)[2]
  cofac <- matrix( 0 , nrow = rows , ncol = cols )
  for (r in range(1,rows)){
    for (c in range(1,cols)){
       e <- (r+c)
       cofac[r,c] = ((-1) ** e) * minor(mat)[r,c]}}
  return(cofac)
}

inverse <- function(mat){
   cofac_t <- t(cofactor(mat))
   inv <- cofac_t /det(mat) 
   return(inv)
}
theta <- c(3,5)
y <- deltasq_Dtheta(c(1,4),2)
y 

inverse(y)
solve(y)

#####################################################################

