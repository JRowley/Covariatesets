Engine <- function(y,d,x,z,beta,xi,gamma,phi,sigma){
  y*d*pmvnorm(lower=c(gamma*x+phi*z,beta*x+xi*d),upper=c(Inf,Inf),sigma=matrix(c(1,sigma,sigma,1),nrow=2,ncol=2)) +
    y*(1-d)*pmvnorm(lower=c(-Inf,beta*x+xi*d),upper=c(gamma*x+phi*z,Inf),sigma=matrix(c(1,sigma,sigma,1),nrow=2,ncol=2)) +
    (1-y)*d*pmvnorm(lower=c(gamma*x+phi*z,-Inf),upper=c(Inf,beta*x+xi*d),sigma=matrix(c(1,sigma,sigma,1),nrow=2,ncol=2)) +
    (1-y)*d*pmvnorm(lower=c(-Inf,-Inf),upper=c(gamma*x+phi*z,beta*x+xi*d),sigma=matrix(c(1,sigma,sigma,1),nrow=2,ncol=2))
}