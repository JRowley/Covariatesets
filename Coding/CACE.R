# ======================================================================== #
# Load data into R
# ======================================================================== #
Source <- url("http://www.ucl.ac.uk/~zctpep9/Data%20archive/AE98Data.RData")
DF <- load(Source)

PUMS80M <- D
rm(list = c("Source","DF","D"))
# ======================================================================== #
# Define RZ
# ======================================================================== #
# Z takes the value 1 if multi2nd == 1
# Z takes the value 0 if multi2nd == 0
PUMS80M$Z <- PUMS80M$multi2nd
# ======================================================================== #
# Define RX
# ======================================================================== #
# X takes the value 1 if whitem == 1
# X takes the value 2 if blackm == 1
# X takes the value 3 if othracem == 1
PUMS80M$X <- PUMS80M$whitem
# ======================================================================== #
# Select relevant data
# ======================================================================== #
Data <- data.frame(Y = PUMS80M$workedm,
                   D = PUMS80M$morekids,
                   X = PUMS80M$X,
                   Z = PUMS80M$Z)
rm(PUMS80M)
# ======================================================================== #
# Define probability space
# ======================================================================== #
library(plyr)
Grid <- expand.grid(unique(Data$Y),
                    unique(Data$D),
                    unique(Data$X),
                    unique(Data$Z))
colnames(Grid) <- c("Y","D","X","Z")
Calc <- function(j){
  W = ifelse(Data$Y == Grid$Y[j] & Data$D == Grid$D[j],1,0)
  V = ifelse(Data$X == Grid$X[j] & Data$Z == Grid$Z[j],1,0)
  R <- lm(W ~ 0 + V)
  return(R$coefficients[[1]])
}
Store <- vector(length = nrow(Grid))
for(j in 1:nrow(Grid)){
  Store[j] = Calc(j)
}
Grid$f <- Store
rm(j)
rm(Store)
# ======================================================================== #
# Define bounds
# ======================================================================== #
Bounds.neg <- expand.grid(unique(Data$X),
                          unique(Data$Z))
colnames(Bounds.neg) <- c("X","Z")
Calc <- function(j){
  A = vector(length = 4)
  # Lower bound Y(0)
  W = ifelse(Data$Y == 0,1,0)
  V = ifelse(Data$X == Bounds.neg$X[j] & Data$Z == Bounds.neg$Z[j],1,0)
  A[1] <- lm(W ~ 0 + V)$coefficients[[1]]
  # Upper bound Y(0)
  W = ifelse(Data$Y == 0 & Data$D == 0,0,1)
  V = ifelse(Data$X == Bounds.neg$X[j] & Data$Z == Bounds.neg$Z[j],1,0)
  A[2] <- lm(W ~ 0 + V)$coefficients[[1]]
  # Lower bound Y(1)
  W = ifelse(Data$Y == 1 & Data$D == 1,1,0)
  V = ifelse(Data$X == Bounds.neg$X[j] & Data$Z == Bounds.neg$Z[j],1,0)
  A[3] <- lm(W ~ 0 + V)$coefficients[[1]]
  # Upper bound Y(1)
  W = ifelse(Data$Y == 1,1,0)
  V = ifelse(Data$X == Bounds.neg$X[j] & Data$Z == Bounds.neg$Z[j],1,0)
  A[4] <- lm(W ~ 0 + V)$coefficients[[1]]
  return(A)
}
Store <- matrix(nrow = nrow(Bounds.neg),ncol = 4)
for(j in 1:nrow(Bounds.neg)){
  Store[j,] = Calc(j)
}
Bounds.neg$l0 <- Store[,1]
Bounds.neg$u0 <- Store[,2]
Bounds.neg$l1 <- Store[,3]
Bounds.neg$u1 <- Store[,4]
rm(j)
rm(Store)
# ======================================================================== #
# Define constraints
# ======================================================================== #
Set <- matrix(nrow = length(unique(Bounds.neg$X)), ncol = 5)
for(j in 1:length(unique(Bounds.neg$X))){
  k = unique(Bounds.neg$X)[j]
  Set[j,1] <- k
  Set[j,2] <- max(Bounds.neg[Bounds.neg$X==k,"l0"])
  Set[j,3] <- min(Bounds.neg[Bounds.neg$X==k,"u0"])
  Set[j,4] <- max(Bounds.neg[Bounds.neg$X==k,"l1"])
  Set[j,5] <- min(Bounds.neg[Bounds.neg$X==k,"u1"])
}
Constraints.neg <- as.data.frame(Set)
rm(Set)
colnames(Constraints.neg) <- c(colnames(Bounds.neg)[-2])
# ======================================================================== #
# Plot
# ======================================================================== #
Plot.data <- Constraints.neg[,1:4]
Plot.data <- melt(Plot.data, id=c("X","l1"))
Plot.data <- Plot.data[,-3]
ggplot(Plot.data, aes(y=l1,x=value,group=X)) + 
  geom_line(aes(col=X))