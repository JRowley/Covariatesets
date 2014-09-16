Calculator <- function(theta){
  beta = theta[1]
  xi = theta[2]
  gamma = theta[3]
  phi = theta[4]
  sig = theta[5]
  require(plyr)
  require(mvtnorm)
  n = nrow(data)
  a = count(data)
  b = count(data,c("X","Z"))
  for(j in 1:16){
    y = data$Y[j]
    d = data$D[j]
    x = data$X[j]
    z = data$Z[j]
    a$p[j] = y*d*(pmvnorm(lower=c(gamma*x+phi*z,beta*x+xi*d),upper=c(Inf,Inf),sigma=matrix(c(1,sig,sig,1),nrow=2,ncol=2))) +
      y*(1-d)*(pmvnorm(lower=c(-Inf,beta*x+xi*d),upper=c(gamma*x+phi*z,Inf),sigma=matrix(c(1,sig,sig,1),nrow=2,ncol=2))) +
      (1-y)*d*(pmvnorm(lower=c(gamma*x+phi*z,-Inf),upper=c(Inf,beta*x+xi*d),sigma=matrix(c(1,sig,sig,1),nrow=2,ncol=2))) +
      (1-y)*(1-d)*(pmvnorm(lower=c(-Inf,-Inf),upper=c(gamma*x+phi*z,beta*x+xi*d),sigma=matrix(c(1,sig,sig,1),nrow=2,ncol=2)))
#     a$o[j] = b[b$X == a$X[j] & b$Z == a$Z[j],"freq"]/n
    a$u[j] = a$p[j]#*a$o[j]
    a$u[j] = log(a$u[j])
    a$out[j] = a$freq[j]*a$u[j]
  }
  find = sum(a$out)
  return(-find)
}

Generator <- function(n,beta,xi,gamma,phi,sig,rho){
  a = rmvnorm(n,mean = c(0,0),sigma=matrix(c(1,rho,rho,1),nrow = 2,ncol = 2))
  a = ifelse(a > 0,1,0)
  b = rmvnorm(n,mean = c(0,0),sigma=matrix(c(1,sig,sig,1),nrow = 2,ncol = 2))
  data = data.frame(X = a[,1],Z = a[,2],U = b[,1],V = b[,2])
  data$D = ifelse(data$X*gamma+data$Z*phi < data$V,1,0)
  data$Y = ifelse(data$X*beta+data$D*gamma < data$U,1,0)
  data = data[,c(1,2,5,6)]
  return(data)
}

n = 10000
beta = -0.5
xi = 1
gamma = 0.1
phi = 1
sig = 0.5
rho = -0.5

data = Generator(n,beta,xi,gamma,phi,sig,rho)

Calculator(c(beta,xi,gamma,phi,sig))

optim(c(beta,xi,gamma,phi,sig),Calculator)
