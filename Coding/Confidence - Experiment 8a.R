load("~/My work/RStudio working directory/Confidence regions for dissertation.RData")
library("mvtnorm", lib.loc="U:/PortableTex/R-3.0.2/library")
library("plyr", lib.loc="U:/PortableTex/R-3.0.2/library")
library("matrixStats", lib.loc="U:/PortableTex/R-3.0.2/library")
rm(Useful.count)
Useful <- PUMS80M[ , c("workedm", "morekids", "multi2nd", "boy1st", "boy2nd", "boys2", "girls2")]
rm(PUMS80M)
# =======================================================================================
# STEP 1
# =======================================================================================
# Set gamma_n.
gamma <- 1 - 0.1/log(nrow(Useful))
# Simulate R draws from the 5-variate standard normal distribution.
R <- 10000000
Draws <- rmvnorm(R, mean = rep(0, 5), sigma = diag(5))
# Set the appropriate quantile of the confidence region.
alpha <- 0.05
# By Bonferroni's inequality, we have that the 'adjusted' level should be 1-alpha/n.
p <- 1-alpha/2
# =======================================================================================
# UPPER BOUND ON RHO.0
# =======================================================================================
# Define a variable Y that takes the value 1 whenever workedm == 0 & morekids == 0.
Useful$Y <- ifelse(Useful$workedm == 0 & Useful$morekids == 0, 1, 0)
# Define variables (X1,X2,X3,X4,X5) that indicate the events (V=1,V=2,V=3,V=4,V=5).
Useful$X1 <- ifelse(Useful$multi2nd == 0 & Useful$boys2 == 1, 1, 0)
Useful$X2 <- ifelse(Useful$multi2nd == 0 & Useful$girls2 == 1, 1, 0)
Useful$X3 <- ifelse(Useful$multi2nd == 0 & Useful$boy1st == 1 & Useful$boy2nd == 0, 1, 0)
Useful$X4 <- ifelse(Useful$multi2nd == 0 & Useful$boy1st == 0 & Useful$boy2nd == 1, 1, 0)
Useful$X5 <- ifelse(Useful$multi2nd == 1, 1, 0)
# Regress Y on (X1,X2,X3) and compute estimated covariance matrix for the parameters.
Regression <- lm(Y ~ 0 + X1 + X2 + X3 + X4 + X5, data = Useful)
# Find variance matrix.
vcov(Regression)
# Adjust the variance matrix by n to get the large sample variance. 
Variance <- vcov(Regression) * nrow(Useful)
# Compute ghat(v) for each v.
g1 <- sqrt(Variance)[, 1]
g2 <- sqrt(Variance)[, 2]
g3 <- sqrt(Variance)[, 3]
g4 <- sqrt(Variance)[, 4]
g5 <- sqrt(Variance)[, 5]
# Compute s(v) for each v.
s1 <- sqrt(sum(g1^2))/sqrt(nrow(Useful))
s2 <- sqrt(sum(g2^2))/sqrt(nrow(Useful))
s3 <- sqrt(sum(g3^2))/sqrt(nrow(Useful))
s4 <- sqrt(sum(g4^2))/sqrt(nrow(Useful))
s5 <- sqrt(sum(g5^2))/sqrt(nrow(Useful))
# Compute the gamma quantile.
g <- data.frame(Draws %*% g1)
colnames(g) <- "v1"
g$v2 <- Draws %*% g2
g$v3 <- Draws %*% g3
g$v4 <- Draws %*% g4
g$v5 <- Draws %*% g5
adjustment <- c(sqrt(sum(g1^2)), sqrt(sum(g2^2)), sqrt(sum(g3^2)), sqrt(sum(g4^2)), sqrt(sum(g5^2)))
# Adjust each column by the norm.
g$v1 <- g$v1/adjustment[1]
g$v2 <- g$v2/adjustment[2]
g$v3 <- g$v3/adjustment[3]
g$v4 <- g$v4/adjustment[4]
g$v5 <- g$v5/adjustment[5]
# Find the maximum over each row.
g <- as.matrix(g)
g <- cbind(g, rowMaxs(g))
# Compute k.auxiliary by finding the gamma quantile of the ordered maxima.
gamma.value <- floor(gamma * R)
k.aux <- sort(g[,3])[gamma.value]
rm(gamma.value)
# Compute the set hat{V}_n.
Check <- data.frame(1-Regression$coefficients)
colnames(Check) <- c("Coefficient")
Check$min.term <- c(
  Check[1,1] + k.aux * s1,
  Check[2,1] + k.aux * s2,
  Check[3,1] + k.aux * s3,
  Check[4,1] + k.aux * s4)
Check$min <- c(
  min(Check$min.term) + 2 * k.aux * s1,
  min(Check$min.term) + 2 * k.aux * s2,
  min(Check$min.term) + 2 * k.aux * s3,
  min(Check$min.term) + 2 * k.aux * s4)
Check$Satisfied <- ifelse(Check$Coefficient <= Check$min, 1, 0)
stop("Code should terminate here and the set hat{V}_n should be found")
# At this point we know which values of v lie in hat{V}_n.
Computation <- matrix(nrow = R, ncol = (Check$Satisfied==1)%*%(rep(1,4)))
j = 1
for(i in 1:nrow(Check)){
  if(Check$Satisfied[i] == 1){
    Computation[,j] <- g[,i]
    j <- j + 1
  }
  else
    print("Not in set")
}
Computation <- cbind(Computation, rowMaxs(Computation))
# Order the final column and find the appropriate quantile.
k.prime <- sort(Computation[,ncol(Computation)])[floor(R*p)]
Bounding <- vector(length = nrow(Check))
for(i in 1:length(Bounding)){
  Bounding[i] <- Check$Coefficient[i] + (1/sqrt(nrow(Useful))) * k.prime * adjustment[i]
}
Bound <- min(Bounding)
U.RHO.0 <- Bound
U.RHO.0
# =======================================================================================
# LOWER BOUND ON RHO.0
# =======================================================================================
load("~/My work/RStudio working directory/Confidence regions for dissertation.RData")
rm(Useful.count)
Useful <- PUMS80M[ , c("workedm", "morekids", "boy1st", "boy2nd", "boys2", "girls2")]
rm(PUMS80M)
# Set gamma_n.
gamma <- 1 - 0.1/log(nrow(Useful))
# Simulate R draws from the 3-variate standard normal distribution.
R <- 10000000
Draws <- rmvnorm(R, mean = rep(0, 4), sigma = diag(4))
# Set the appropriate quantile of the confidence region.
alpha <- 0.05
# By Bonferroni's inequality, we have that the 'adjusted' level should be 1-alpha/n.
p <- 1-alpha/2
# Define a variable Y that takes the value 1 whenever workedm == 1.
Useful$Y <- ifelse(Useful$workedm == 1, 1, 0)
# Define variables (X1,X2,X3,X4,X5) that indicate the events (V=1,V=2,V=3,V=4,V=5).
Useful$X1 <- ifelse(Useful$multi2nd == 0 & Useful$boys2 == 1, 1, 0)
Useful$X2 <- ifelse(Useful$multi2nd == 0 & Useful$girls2 == 1, 1, 0)
Useful$X3 <- ifelse(Useful$multi2nd == 0 & Useful$boy1st == 1 & Useful$boy2nd == 0, 1, 0)
Useful$X4 <- ifelse(Useful$multi2nd == 0 & Useful$boy1st == 0 & Useful$boy2nd == 1, 1, 0)
Useful$X5 <- ifelse(Useful$multi2nd == 1, 1, 0)
# Regress Y on (X1,X2,X3) and compute estimated covariance matrix for the parameters.
Regression <- lm(Y ~ 0 + X1 + X2 + X3 + X4 + X5, data = Useful)
# Find variance matrix.
vcov(Regression)
# Adjust the variance matrix by n to get the large sample variance. 
Variance <- vcov(Regression) * nrow(Useful)
# Compute ghat(v) for each v.
g1 <- sqrt(Variance)[, 1]
g2 <- sqrt(Variance)[, 2]
g3 <- sqrt(Variance)[, 3]
g4 <- sqrt(Variance)[, 4]
g5 <- sqrt(Variance)[, 5]
# Compute s(v) for each v.
s1 <- sqrt(sum(g1^2))/sqrt(nrow(Useful))
s2 <- sqrt(sum(g2^2))/sqrt(nrow(Useful))
s3 <- sqrt(sum(g3^2))/sqrt(nrow(Useful))
s4 <- sqrt(sum(g4^2))/sqrt(nrow(Useful))
s5 <- sqrt(sum(g5^2))/sqrt(nrow(Useful))
# Compute the gamma quantile.
g <- data.frame(Draws %*% g1)
colnames(g) <- "v1"
g$v2 <- Draws %*% g2
g$v3 <- Draws %*% g3
g$v4 <- Draws %*% g4
g$v5 <- Draws %*% g5
adjustment <- c(sqrt(sum(g1^2)), sqrt(sum(g2^2)), sqrt(sum(g3^2)), sqrt(sum(g4^2)), sqrt(sum(g5^2)))
# Adjust each column by the norm.
g$v1 <- g$v1/adjustment[1]
g$v2 <- g$v2/adjustment[2]
g$v3 <- g$v3/adjustment[3]
g$v4 <- g$v4/adjustment[4]
g$v5 <- g$v5/adjustment[5]
# Find the maximum over each row.
g <- as.matrix(g)
g <- cbind(g, rowMaxs(g))
# Compute k.auxiliary by finding the gamma quantile of the ordered maxima.
gamma.value <- floor(gamma * R)
k.aux <- sort(g[,3])[gamma.value]
rm(gamma.value)
# Compute the set hat{V}_n.
Check <- data.frame(-Regression$coefficients)
colnames(Check) <- c("Coefficient")
Check$min.term <- c(
  Check[1,1] + k.aux * s1,
  Check[2,1] + k.aux * s2,
  Check[3,1] + k.aux * s3,
  Check[4,1] + k.aux * s4)
Check$min <- c(
  min(Check$min.term) + 2 * k.aux * s1,
  min(Check$min.term) + 2 * k.aux * s2,
  min(Check$min.term) + 2 * k.aux * s3,
  min(Check$min.term) + 2 * k.aux * s4)
Check$Satisfied <- ifelse(Check$Coefficient <= Check$min, 1, 0)
stop("Code should terminate here and the set hat{V}_n should be found")
# At this point we know which values of v lie in hat{V}_n.
Computation <- matrix(nrow = R, ncol = (Check$Satisfied==1)%*%(rep(1,4)))
j = 1
for(i in 1:nrow(Check)){
  if(Check$Satisfied[i] == 1){
    Computation[,j] <- g[,i]
    j <- j + 1
  }
  else
    print("Not in set")
}
Computation <- cbind(Computation, rowMaxs(Computation))
# Order the final column and find the appropriate quantile.
k.prime <- sort(Computation[,ncol(Computation)])[floor(R*p)]
Bounding <- vector(length = nrow(Check))
for(i in 1:length(Bounding)){
  Bounding[i] <- Check$Coefficient[i] + (1/sqrt(nrow(Useful))) * k.prime * adjustment[i] 
}
Bound <- min(Bounding)
L.RHO.0 <- -1 * Bound
L.RHO.0
# =======================================================================================
# UPPER BOUND ON RHO.0+RHO.1
# =======================================================================================
load("~/My work/RStudio working directory/Confidence regions for dissertation.RData")
rm(Useful.count)
Useful <- PUMS80M[ , c("workedm", "morekids", "boy1st", "boy2nd", "boys2", "girls2")]
rm(PUMS80M)
# Set gamma_n.
gamma <- 1 - 0.1/log(nrow(Useful))
# Simulate R draws from the 3-variate standard normal distribution.
R <- 10000000
Draws <- rmvnorm(R, mean = rep(0, 4), sigma = diag(4))
# Set the appropriate quantile of the confidence region.
alpha <- 0.05
# By Bonferroni's inequality, we have that the 'adjusted' level should be 1-alpha/n.
p <- 1-alpha/2
# Define a variable Y that takes the value 1 whenever workedm == 1.
Useful$Y <- ifelse(Useful$workedm == 1, 1, 0)
# Define variables (X1,X2,X3,X4,X5) that indicate the events (V=1,V=2,V=3,V=4,V=5).
Useful$X1 <- ifelse(Useful$multi2nd == 0 & Useful$boys2 == 1, 1, 0)
Useful$X2 <- ifelse(Useful$multi2nd == 0 & Useful$girls2 == 1, 1, 0)
Useful$X3 <- ifelse(Useful$multi2nd == 0 & Useful$boy1st == 1 & Useful$boy2nd == 0, 1, 0)
Useful$X4 <- ifelse(Useful$multi2nd == 0 & Useful$boy1st == 0 & Useful$boy2nd == 1, 1, 0)
Useful$X5 <- ifelse(Useful$multi2nd == 1, 1, 0)
# Regress Y on (X1,X2,X3) and compute estimated covariance matrix for the parameters.
Regression <- lm(Y ~ 0 + X1 + X2 + X3 + X4 + X5, data = Useful)
# Find variance matrix.
vcov(Regression)
# Adjust the variance matrix by n to get the large sample variance. 
Variance <- vcov(Regression) * nrow(Useful)
# Compute ghat(v) for each v.
g1 <- sqrt(Variance)[, 1]
g2 <- sqrt(Variance)[, 2]
g3 <- sqrt(Variance)[, 3]
g4 <- sqrt(Variance)[, 4]
g5 <- sqrt(Variance)[, 5]
# Compute s(v) for each v.
s1 <- sqrt(sum(g1^2))/sqrt(nrow(Useful))
s2 <- sqrt(sum(g2^2))/sqrt(nrow(Useful))
s3 <- sqrt(sum(g3^2))/sqrt(nrow(Useful))
s4 <- sqrt(sum(g4^2))/sqrt(nrow(Useful))
s5 <- sqrt(sum(g5^2))/sqrt(nrow(Useful))
# Compute the gamma quantile.
g <- data.frame(Draws %*% g1)
colnames(g) <- "v1"
g$v2 <- Draws %*% g2
g$v3 <- Draws %*% g3
g$v4 <- Draws %*% g4
g$v5 <- Draws %*% g5
adjustment <- c(sqrt(sum(g1^2)), sqrt(sum(g2^2)), sqrt(sum(g3^2)), sqrt(sum(g4^2)), sqrt(sum(g5^2)))
# Adjust each column by the norm.
g$v1 <- g$v1/adjustment[1]
g$v2 <- g$v2/adjustment[2]
g$v3 <- g$v3/adjustment[3]
g$v4 <- g$v4/adjustment[4]
g$v5 <- g$v5/adjustment[5]
# Find the maximum over each row.
g <- as.matrix(g)
g <- cbind(g, rowMaxs(g))
# Compute k.auxiliary by finding the gamma quantile of the ordered maxima.
gamma.value <- floor(gamma * R)
k.aux <- sort(g[,3])[gamma.value]
rm(gamma.value)
# Compute the set hat{V}_n.
Check <- data.frame(Regression$coefficients)
colnames(Check) <- c("Coefficient")
Check$min.term <- c(
  Check[1,1] + k.aux * s1,
  Check[2,1] + k.aux * s2,
  Check[3,1] + k.aux * s3,
  Check[4,1] + k.aux * s4)
Check$min <- c(
  min(Check$min.term) + 2 * k.aux * s1,
  min(Check$min.term) + 2 * k.aux * s2,
  min(Check$min.term) + 2 * k.aux * s3,
  min(Check$min.term) + 2 * k.aux * s4)
Check$Satisfied <- ifelse(Check$Coefficient <= Check$min, 1, 0)
stop("Code should terminate here and the set hat{V}_n should be found")
# At this point we know which values of v lie in hat{V}_n.
Computation <- matrix(nrow = R, ncol = (Check$Satisfied==1)%*%(rep(1,4)))
j = 1
for(i in 1:nrow(Check)){
  if(Check$Satisfied[i] == 1){
    Computation[,j] <- g[,i]
    j <- j + 1
  }
  else
    print("Not in set")
}
Computation <- cbind(Computation, rowMaxs(Computation))
# Order the final column and find the appropriate quantile.
k.prime <- sort(Computation[,ncol(Computation)])[floor(R*p)]
Bounding <- vector(length = nrow(Check))
for(i in 1:length(Bounding)){
  Bounding[i] <- Check$Coefficient[i] + (1/sqrt(nrow(Useful))) * k.prime * adjustment[i]
}
Bound <- min(Bounding)
U.RHO.01 <- Bound
U.RHO.01
# =======================================================================================
# LOWER BOUND ON RHO.0+RHO.1
# =======================================================================================
load("~/My work/RStudio working directory/Confidence regions for dissertation.RData")
rm(Useful.count)
Useful <- PUMS80M[ , c("workedm", "morekids", "boy1st", "boy2nd", "boys2", "girls2")]
rm(PUMS80M)
# Set gamma_n.
gamma <- 1 - 0.1/log(nrow(Useful))
# Simulate R draws from the 3-variate standard normal distribution.
R <- 10000000
Draws <- rmvnorm(R, mean = rep(0, 4), sigma = diag(4))
# Set the appropriate quantile of the confidence region.
alpha <- 0.05
# By Bonferroni's inequality, we have that the 'adjusted' level should be 1-alpha/n.
p <- 1-alpha/2
# Define a variable Y that takes the value 1 whenever workedm == 1 & morekids == 1.
Useful$Y <- ifelse(Useful$workedm == 1 & Useful$morekids == 1, 1, 0)
# Define variables (X1,X2,X3,X4,X5) that indicate the events (V=1,V=2,V=3,V=4,V=5).
Useful$X1 <- ifelse(Useful$multi2nd == 0 & Useful$boys2 == 1, 1, 0)
Useful$X2 <- ifelse(Useful$multi2nd == 0 & Useful$girls2 == 1, 1, 0)
Useful$X3 <- ifelse(Useful$multi2nd == 0 & Useful$boy1st == 1 & Useful$boy2nd == 0, 1, 0)
Useful$X4 <- ifelse(Useful$multi2nd == 0 & Useful$boy1st == 0 & Useful$boy2nd == 1, 1, 0)
Useful$X5 <- ifelse(Useful$multi2nd == 1, 1, 0)
# Regress Y on (X1,X2,X3) and compute estimated covariance matrix for the parameters.
Regression <- lm(Y ~ 0 + X1 + X2 + X3 + X4 + X5, data = Useful)
# Find variance matrix.
vcov(Regression)
# Adjust the variance matrix by n to get the large sample variance. 
Variance <- vcov(Regression) * nrow(Useful)
# Compute ghat(v) for each v.
g1 <- sqrt(Variance)[, 1]
g2 <- sqrt(Variance)[, 2]
g3 <- sqrt(Variance)[, 3]
g4 <- sqrt(Variance)[, 4]
g5 <- sqrt(Variance)[, 5]
# Compute s(v) for each v.
s1 <- sqrt(sum(g1^2))/sqrt(nrow(Useful))
s2 <- sqrt(sum(g2^2))/sqrt(nrow(Useful))
s3 <- sqrt(sum(g3^2))/sqrt(nrow(Useful))
s4 <- sqrt(sum(g4^2))/sqrt(nrow(Useful))
s5 <- sqrt(sum(g5^2))/sqrt(nrow(Useful))
# Compute the gamma quantile.
g <- data.frame(Draws %*% g1)
colnames(g) <- "v1"
g$v2 <- Draws %*% g2
g$v3 <- Draws %*% g3
g$v4 <- Draws %*% g4
g$v5 <- Draws %*% g5
adjustment <- c(sqrt(sum(g1^2)), sqrt(sum(g2^2)), sqrt(sum(g3^2)), sqrt(sum(g4^2)), sqrt(sum(g5^2)))
# Adjust each column by the norm.
g$v1 <- g$v1/adjustment[1]
g$v2 <- g$v2/adjustment[2]
g$v3 <- g$v3/adjustment[3]
g$v4 <- g$v4/adjustment[4]
g$v5 <- g$v5/adjustment[5]
# Find the maximum over each row.
g <- as.matrix(g)
g <- cbind(g, rowMaxs(g))
# Compute k.auxiliary by finding the gamma quantile of the ordered maxima.
gamma.value <- floor(gamma * R)
k.aux <- sort(g[,3])[gamma.value]
rm(gamma.value)
# Compute the set hat{V}_n.
Check <- data.frame(-Regression$coefficients)
colnames(Check) <- c("Coefficient")
Check$min.term <- c(
  Check[1,1] + k.aux * s1,
  Check[2,1] + k.aux * s2,
  Check[3,1] + k.aux * s3,
  Check[3,1] + k.aux * s3)
Check$min <- c(
  min(Check$min.term) + 2 * k.aux * s1,
  min(Check$min.term) + 2 * k.aux * s2,
  min(Check$min.term) + 2 * k.aux * s3,
  min(Check$min.term) + 2 * k.aux * s4)
Check$Satisfied <- ifelse(Check$Coefficient <= Check$min, 1, 0)
stop("Code should terminate here and the set hat{V}_n should be found")
# At this point we know which values of v lie in hat{V}_n.
Computation <- matrix(nrow = R, ncol = (Check$Satisfied==1)%*%(rep(1,4)))
j = 1
for(i in 1:nrow(Check)){ 
  if(Check$Satisfied[i] == 1){
    Computation[,j] <- g[,i]
    j <- j + 1
  }
  else
    print("Not in set")
}
Computation <- cbind(Computation, rowMaxs(Computation))
# Order the final column and find the appropriate quantile.
k.prime <- sort(Computation[,ncol(Computation)])[floor(R*p)]
Bounding <- vector(length = nrow(Check))
for(i in 1:length(Bounding)){
  Bounding[i] <- Check$Coefficient[i] + (1/sqrt(nrow(Useful))) * k.prime * adjustment[i]
}
Bound <- min(Bounding)
L.RHO.01 <- -1 * Bound
L.RHO.01
# =======================================================================================
# RESULTS
# =======================================================================================
rho.0_lowerbound
rho.0_upperbound
L.RHO.0
U.RHO.0
rho.0andrho.1_lowerbound
rho.0andrho.1_upperbound
L.RHO.01
U.RHO.01