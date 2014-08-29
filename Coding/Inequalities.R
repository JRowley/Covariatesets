ptm <- proc.time()
# ======================================================================= #
# Setup.
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

Source <- url("http://www.ucl.ac.uk/~zctpep9/Data%20archive/AE98Data.RData")
DF <- load(Source)

PUMS80M <- D
rm(list = c("Source","DF","D"))

Data <- PUMS80M[,c("workedm","morekids","hispm","multi2nd")]
colnames(Data) <- c("Y","D","X","Z")

rm(PUMS80M)

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
rm(Data)

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
    r[[k]] <- ifelse(prod(s$upper >= s$lower) == 0,NA,r[k])
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

bound <- a(bound)

rm(a)

proc.time() - ptm