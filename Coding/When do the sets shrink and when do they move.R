ptm <- proc.time()

library(mvtnorm)
library(iterpc)

Generator <- function(n,d,x.y,x.d,z,sigma,mu.x,mu.z,co){
  Z <- c(rep(c(0,1),0.5*n))
  X <- c(Z[1:(n*0.5*(1+co))],1-Z[(n*0.5*(1+co)+1):n])
  d
  set.seed(10)
  e <- rmvnorm(n,mean=c(0,0.5*mu.x+0.5*mu.z),sigma=matrix(c(1,sigma,sigma,1),nrow=2,ncol=2))
  D <- ifelse(x.d*X+z*Z<e[,2],1,0)
  Y <- ifelse(x.y*X+d*D<e[,1],1,0)
  data <- data.frame(workedm=Y,morekids=D,hispm=X,multi2nd=Z)
  return(data)
}

# n = 1e6
# d = 1
# x.y = -0.5
# x.d = 0.75
# z = -1
# sigma = 0.5
# mu.x = 1
# mu.z = -0.5
# co = 0.3

n = 1e3
d = 1
x.y = -0.5
x.d = 0.75
z = -1
sigma = 0.5
mu.x = 0
mu.z = 0
co = 0

PUMS80M <- Generator(n,d,x.y,x.d,z,sigma,mu.x,mu.z,co)
# ======================================================================= #
# Setup.
# ======================================================================= #

tol = 1e-8

D <- PUMS80M

library(iterpc)

# Get all possible strict orderings of r(d).
a = iterpc(2, 2, ordered = T, replace = T)
Order <- data.frame(getall(a))
Output <- matrix(nrow=nrow(Order),ncol=2)
for(j in 1:nrow(Order)){
  Output[j,] = as.numeric(factor(as.numeric(Order[j,])))
}
Order <- unique(Output)
rm(Output)
rm(a)
# Rename columns of Order.
a = 0:1 
a = as.vector(paste("f",a,sep="."))
colnames(Order) <- a
rm(a)
rm(j)

# ======================================================================= #
# Get a grid for (Y,D,X,Z) and fill in probabilities.
# ======================================================================= #

Data <- D[,c("workedm","morekids","multi2nd")]
colnames(Data) <- c("Y","D","Z")

rm(PUMS80M)

P <- expand.grid(Y=unique(Data$Y),
                 D=unique(Data$D),
                 Z=unique(Data$Z))

library(plyr)

a <- count(Data)
b <- count(Data,c("Z"))
colnames(a) <- c(colnames(a[,1:3]),"Total")

Prob <- merge(a,b,by=c("Z"),all=T)
Prob$p <- Prob$Total/Prob$freq
Prob <- Prob[,c(1:3,6)]
rm(a)
rm(b)

P <- merge(P,Prob,by=c("Y","D","Z"),all=T)
P[is.na(P)] = 0
rm(Prob)
rm(Data)

# ====================================================================== #
# Compute the lower probability for an order.
# ====================================================================== #

grid <- expand.grid(d = unique(P$D),
                    z = unique(P$Z))

grid$reference <- 1 + grid$d 

# Let j be the ordering that is being considered.
# Let k be the row of grid that is being considered.

# I will take a particular ordering and then create a list for that 
# ordering. The list will consist of one slice for each row of the grid,
# and will be the set A that I have defined in the text.


a = function(j){
  output <- list()
  for(k in 1:nrow(grid)){
    output[[k]] <- which(Order[j,] <=
                           Order[j,grid$reference[k]]) - 1        
  }
  return(output)
}

b = function(j){
  intermediate <- a(j)
  output <- vector(length=nrow(grid))
  for(k in 1:nrow(grid)){
    output[k] <- sum(P[P$Y == 0 & P$D %in% unlist(intermediate[k]) &
                         P$Z == grid$z[k],]$p)    
  }
  return(output)
}

lower = function(){
  output <- list()
  for(j in 1:nrow(Order)){
    intermediate <- cbind(grid,b(j))
    colnames(intermediate) <- c(colnames(grid),"inequal")
    intermediate <- ddply(intermediate,.(d),summarize,bound=max(inequal))
    output[[j]] <- intermediate
  }
  return(output)
}

bound.l <- lower()

rm(list=c("a","b","lower"))

# ====================================================================== #
# Compute the upper probability for an order.
# ====================================================================== #

a = function(j){
  output <- list()
  for(k in 1:nrow(grid)){
    output[[k]] <- which(Order[j,] >=
                           Order[j,grid$reference[k]]) - 1        
  }
  return(output)
}

b = function(j){
  intermediate <- a(j)
  output <- vector(length=nrow(grid))
  for(k in 1:nrow(grid)){
    output[k] <- sum(P[P$Y == 1 & P$D %in% unlist(intermediate[k]) &
                         P$Z == grid$z[k],]$p)    
  }
  return(output)
}

upper = function(){
  output <- list()
  for(j in 1:nrow(Order)){
    intermediate <- cbind(grid,b(j))
    colnames(intermediate) <- c(colnames(grid),"inequal")
    intermediate <- ddply(intermediate,.(d),summarize,bound=1-max(inequal))
    output[[j]] <- intermediate
  }
  return(output)
}

bound.u <- upper()

rm(list=c("a","b","upper"))

# ====================================================================== #
# Store them all in one list.
# ====================================================================== #

a = function(r,s){
  bound <- list()
  for(k in 1:length(r)){
    bound[[k]] <- merge(r[k],s[k],by=c("d"),all=T)
    colnames(bound[[k]]) <- c("d","lower","upper")
  }
  return(bound)
}

b = function(r){
  for(k in 1:length(r)){
    s = data.frame(r[k])
    r[[k]] <- ifelse(prod(s$upper >= (s$lower-tol)) == 0,NA,r[k])
  }
  return(r)
}

bound <- a(bound.l,bound.u)
bound <- b(bound)

rm(list=c("a","b"))

# ====================================================================== #
# Remove bounds that do not satisfy the ordering.
# ====================================================================== #

a = function(r){
  Store <- list()
  for(j in 1:length(r[is.na(r)==F])){
    s <- data.frame(r[is.na(r)==F][j])
    s <- s[with(s, order(d)),]
    s <- cbind(s,Order[is.na(r)==F,])
    colnames(s) <- c("d","lower","upper","order")
    Store[[j]] <- s
  }
  return(Store)
}

bound.x1 <- a(bound)

rm(a)

# ======================================================================= #
# ======================================================================= #
# ======================================================================= #
# Do again for X=2.
# ======================================================================= #
# ======================================================================= #
# ======================================================================= #

# Specify the support of X.
X = 2
# I will consider whether a mother is white or not.

library(iterpc)

# Get all possible strict orderings of p(d,x).
a = iterpc(2*X, 2*X, ordered = T, replace = T)
Order <- data.frame(getall(a))
Output <- matrix(nrow=nrow(Order),ncol=4)
for(j in 1:nrow(Order)){
  Output[j,] = as.numeric(factor(as.numeric(Order[j,])))
}
Order <- unique(Output)
rm(Output)
rm(a)
# Rename columns of Order.
a = as.vector(outer(0:1, 1:X-1, paste, sep="")) 
a = as.vector(paste("f",a,sep="."))
colnames(Order) <- a
rm(a)
# Redefine X in letters to match naming convention of Order.
X = letters[1:X]

# ======================================================================= #
# Get a grid for (Y,D,X,Z) and fill in probabilities.
# ======================================================================= #

Data <- D[,c("workedm","morekids","hispm","multi2nd")]
colnames(Data) <- c("Y","D","X","Z")

P <- expand.grid(Y=unique(Data$Y),
                 D=unique(Data$D),
                 X=unique(Data$X),
                 Z=unique(Data$Z))

library(plyr)

a <- count(Data)
b <- count(Data,c("X","Z"))
colnames(a) <- c(colnames(a[,1:4]),"Total")

Prob <- merge(a,b,by=c("X","Z"),all=T)
Prob$p <- Prob$Total/Prob$freq
Prob <- Prob[,c(1:4,7)]
rm(a)
rm(b)

P <- merge(P,Prob,by=c("Y","D","X","Z"),all=T)
P[is.na(P)] = 0
rm(Prob)

# ====================================================================== #
# Compute the lower probability for an order.
# ====================================================================== #

grid <- expand.grid(d = unique(P$D),
                    x = unique(P$X),
                    eta = unique(P$X),
                    z = unique(P$Z))

grid$reference <- 1 + grid$d + 2*grid$x

# Let j be the ordering that is being considered.
# Let k be the row of grid that is being considered.

# I will take a particular ordering and then create a list for that 
# ordering. The list will consist of one slice for each row of the grid,
# and will be the set A that I have defined in the text.


a = function(j){
  output <- list()
  for(k in 1:nrow(grid)){
    q = grid[k,]$eta
    output[[k]] <- which(Order[j,((1+2*q):(2+2*q))] <=
                           Order[j,grid$reference[k]]) - 1        
  }
  return(output)
}

b = function(j){
  intermediate <- a(j)
  output <- vector(length=nrow(grid))
  for(k in 1:nrow(grid)){
    output[k] <- sum(P[P$Y == 0 & P$D %in% unlist(intermediate[k]) &
                         P$X == grid$eta[k] & P$Z == grid$z[k],]$p)    
  }
  return(output)
}

lower = function(){
  output <- list()
  for(j in 1:nrow(Order)){
    intermediate <- cbind(grid,b(j))
    colnames(intermediate) <- c(colnames(grid),"inequal")
    intermediate <- ddply(intermediate,.(d,x),summarize,bound=max(inequal))
    output[[j]] <- intermediate
  }
  return(output)
}

bound.l <- lower()

rm(list=c("a","b","lower"))

# ====================================================================== #
# Compute the upper probability for an order.
# ====================================================================== #

a = function(j){
  output <- list()
  for(k in 1:nrow(grid)){
    q = grid[k,]$eta
    output[[k]] <- which(Order[j,((1+2*q):(2+2*q))] >=
                           Order[j,grid$reference[k]]) - 1        
  }
  return(output)
}

b = function(j){
  intermediate <- a(j)
  output <- vector(length=nrow(grid))
  for(k in 1:nrow(grid)){
    output[k] <- sum(P[P$Y == 1 & P$D %in% unlist(intermediate[k]) &
                         P$X == grid$eta[k] & P$Z == grid$z[k],]$p)    
  }
  return(output)
}

upper = function(){
  output <- list()
  for(j in 1:nrow(Order)){
    intermediate <- cbind(grid,b(j))
    colnames(intermediate) <- c(colnames(grid),"inequal")
    intermediate <- ddply(intermediate,.(d,x),summarize,bound=1-max(inequal))
    output[[j]] <- intermediate
  }
  return(output)
}

bound.u <- upper()

rm(list=c("a","b","upper"))

# ====================================================================== #
# Store them all in one list.
# ====================================================================== #

a = function(r,s){
  bound <- list()
  for(k in 1:length(r)){
    bound[[k]] <- merge(r[k],s[k],by=c("d","x"),all=T)
    colnames(bound[[k]]) <- c("d","x","lower","upper")
  }
  return(bound)
}

b = function(r){
  for(k in 1:length(r)){
    s = data.frame(r[k])
    r[[k]] <- ifelse(prod(s$upper >= (s$lower-tol)) == 0,NA,r[k])
  }
  return(r)
}

bound <- a(bound.l,bound.u)
bound <- b(bound)

rm(list=c("a","b"))

# ====================================================================== #
# Remove bounds that do not satisfy the ordering.
# ====================================================================== #

a = function(r){
  Store <- list()
  for(j in 1:length(r[is.na(r)==F])){
    s <- data.frame(r[is.na(r)==F][j])
    s <- s[with(s, order(x,d)),]
    s <- cbind(s,Order[is.na(r)==F,][j,])
    colnames(s) <- c("d","x","lower","upper","order")
    Store[[j]] <- s
  }
  return(Store)
}

b = function(r){
  Store <- list()
  for(j in 1:length(r[is.na(r)==F])){
    s <- data.frame(r[is.na(r)==F][j])
    s <- s[with(s, order(x,d)),]
    s <- cbind(s,Order[is.na(r)==F,])
    colnames(s) <- c("d","x","lower","upper","order")
    Store[[j]] <- s
  }
  return(Store)
}

bound <- if(length(bound[is.na(bound)==F])>1){
  a(bound)} else{
    b(bound)}

rm(a)
rm(b)

# ====================================================================== #
# Remove bounds that do not satisfy the ordering (externally invalid).
# ====================================================================== #

# a = function(r){
#   Store <- vector(length = length(r[is.na(r)==F]))
#   for(j in 1:length(r[is.na(r)==F])){
#     s <- data.frame(r[is.na(r)==F][j])
#     s <- arrange(s,order)
#     bash <- vector(length = nrow(s))
#     for(k in 1:nrow(s)){
#       bash[k] <- prod(s$upper[k]<=s[s$order>s$order[k],3])
#     }
#     Store[j] = prod(bash)
#   }
#   return(Store)
# }
# 
# b <- which(a(bound)==1)
# bound <- bound[b]
# 
# rm(a)
# rm(b)

# ====================================================================== #
# Find the weighted bound.
# ====================================================================== #

a <- count(Data,"X")
a$p <- a$freq/sum(a$freq)
a <- a[,c(1,3)]
colnames(a) <- c("x","p")

b = function(r,s){
  Store <- list()
  for(j in 1:length(r)){
    q = merge(data.frame(r[j]),s)
    q = cbind(q[,1:4],q[,6])
    colnames(q) <- c("x","d","l","u","p")
    q$lower <- q$p*q$l
    q$upper <- q$p*q$u
    q <- q[,c(1,2,6,7)]
    q <- ddply(q,.(d),summarize,lower=sum(lower),upper=sum(upper))
    Store[[j]] <- q
  }
  return(Store)
}

Bound <- b(bound,a)

rm(a)
rm(b)

# ====================================================================== #
# Compute ACE(D->Y|x).
# ====================================================================== #

a = function(r){
  Store <- list()
  for(j in 1:length(r)){
    q = data.frame(r[j])
    b <- q[q$d==0 & q$x==0,]$lower-q[q$d==1 & q$x==0,]$upper
    f <- q[q$d==0 & q$x==0,]$upper-q[q$d==1 & q$x==0,]$lower
    g <- q[q$d==0 & q$x==1,]$lower-q[q$d==1 & q$x==1,]$upper
    h <- q[q$d==0 & q$x==1,]$upper-q[q$d==1 & q$x==1,]$lower
    Store[[j]] <- list(c(b,f),c(g,h))
  }
  return(Store)
}

a(bound)

rm(a)

# ====================================================================== #
# Compute ACE(D->Y).
# ====================================================================== #

a = function(r){
  Store <- list()
  for(j in 1:length(r)){
    q = data.frame(r[j])
    b <- q[q$d==0,]$lower-q[q$d==1,]$upper
    f <- q[q$d==0,]$upper-q[q$d==1,]$lower
    Store[[j]] <- c(b,f)
  }
  return(Store)
}

a(Bound)
a(bound.x1)
# bound
# bound.x1
store <- as.numeric(a(bound.x1)[[1]])
Store <- matrix(unlist(a(Bound)), ncol = 2, byrow = TRUE)

rm(a)

c(min(store) - min(Store),max(store) - max(Store))

proc.time() - ptm