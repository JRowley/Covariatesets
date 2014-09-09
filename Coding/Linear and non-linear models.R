n = 200000

Data <- function(n,end,cor,aff){
  A <- data.frame(rmvnorm(n,mean=c(0,0),sigma=matrix(c(1,end,end,1),nrow=2,ncol=2)))
  colnames(A) <- c("U","V")
  B <- data.frame(rmvnorm(n,mean=c(0,0),sigma=matrix(c(1,cor,cor,1),nrow=2,ncol=2)))
  colnames(B) <- c("X","Z")
  D <- data.frame(A,B)
  D$X <- ifelse(D$X>0,1,0)
  D$Z <- ifelse(D$Z>0,1,0)
  D$d <- ifelse(D$V>0.5*D$Z+aff*0.5*D$X,1,0)
  D$y <- ifelse(D$U>0.5*D$d+0.5*D$X,1,0)
  return(D)
}

#Spec1 - No endog, X does not affect D, zero corr X and Z
Frame <- Data(n,0,0,0)
tsls(y~d+X,~X+Z,data=Frame)
tsls(y~d,~Z,data=Frame)
# Doesn't matter whether I put in X or not

#Spec2 - endog, X does not affect D, zero corr X and Z
Frame <- Data(n,0.5,0,0)
tsls(y~d+X,~X+Z,data=Frame)
tsls(y~d,~Z,data=Frame)
# Again - does not matter

#Spec3 - No endog, X does affect D, zero corr X and Z
Frame <- Data(n,0,0,1)
tsls(y~d+X,~X+Z,data=Frame)
tsls(y~d,~Z,data=Frame)
# Again - does not matter

#Spec4 - No endog, X does not affect D, corr X and Z
Frame <- Data(n,0,0.5,0)
tsls(y~d+X,~X+Z,data=Frame)
tsls(y~d,~Z,data=Frame)
# Does matter

#Spec5 - endog, X does affect D, zero corr X and Z
Frame <- Data(n,0.5,0,1)
tsls(y~d+X,~X+Z,data=Frame)
tsls(y~d,~Z,data=Frame)
# Does not matter

#Spec6 - endog, X does not affect D, corr X and Z
Frame <- Data(n,0.5,0.5,0)
tsls(y~d+X,~X+Z,data=Frame)
tsls(y~d,~Z,data=Frame)
# Does matter

#Spec7 - endog, X does affect D, corr X and Z
Frame <- Data(n,0.5,0.5,1)
tsls(y~d+X,~X+Z,data=Frame)
tsls(y~d,~Z,data=Frame)
# Does matter

#Spec8 - no endog, X does affect D, corr X and Z
Frame <- Data(n,0,0.5,1)
tsls(y~d+X,~X+Z,data=Frame)
tsls(y~d,~Z,data=Frame)
# Does matter