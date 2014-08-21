Utype <- 1:16
Vtype <- 1:16
Xtype <- 0:1
Ztype <- 0:1
Unobservable <- expand.grid(Utype,Vtype)
Observable <- expand.grid(Xtype,Ztype)
colnames(Unobservable) <- c("Utype","Vtype")
colnames(Observable) <- c("X","Z")
rm(list=c("Utype","Vtype","Xtype","Ztype"))
# =================================================================================================== #
# Probability of events  
# =================================================================================================== #
Unobservable$p <- runif(nrow(Unobservable)) 
Unobservable$p <- Unobservable$p/sum(Unobservable$p)
Observable$p <- runif(nrow(Observable))
Observable$p <- Observable$p/sum(Observable$p)
# =================================================================================================== #
# Create new data frame 
# =================================================================================================== #
Event <- expand.grid(unique(Unobservable$Utype),
                     unique(Unobservable$Vtype),
                     unique(Observable$X),
                     unique(Observable$Z))
colnames(Event) <- c("Utype","Vtype","X","Z")
# =================================================================================================== #
# Get probabilities
# =================================================================================================== #
Calculator <- function(j){
  a = Unobservable$p[Unobservable$Utype == Event$Utype[j] & Unobservable$Vtype == Event$Vtype[j]] 
  b = Observable$p[Observable$X == Event$X[j] & Observable$Z == Event$Z[j]]
  return(a*b)
}
Event$p <- rep(0,nrow(Event))
for(j in 1:nrow(Event)){
  Event$p[j] <- Calculator(j)
}
rm(j)
rm(Calculator)
# =================================================================================================== #
# Description of typespace
# =================================================================================================== #
# Utype

# Type 1 : Never taker, x=0
# Type 2 : Never taker, x=1
# Type 3 : Complier, x=0
# Type 4 : Complier, x=1
# Type 5 : Defier, x=0
# Type 6 : Defier, x=1
# Type 7 : Always taker, x=0
# Type 8 : Always taker, x=1

# Vtype

# Type 1 : Never taker, x=0
# Type 2 : Never taker, x=1
# Type 3 : Complier, x=0
# Type 4 : Complier, x=1
# Type 5 : Defier, x=0
# Type 6 : Defier, x=1
# Type 7 : Always taker, x=0
# Type 8 : Always taker, x=1
# =================================================================================================== #
# Y and D
# =================================================================================================== #
Calculator <- function(j){
  a = (1-Event$X[j])*(ifelse(Event$Vtype[j] %in% 1:4,0,
                          ifelse(Event$Vtype[j] %in% 5:8,Event$Z[j],
                                 ifelse(Event$Vtype[j] %in% 9:12,1-Event$Z[j],1)))) +
    (Event$X[j])*(ifelse(Event$Vtype[j] %in% c(1,5,9,13),0,
                         ifelse(Event$Vtype[j] %in% c(2,6,10,14),Event$Z[j],
                                ifelse(Event$Vtype[j] %in% c(3,7,11,15),1-Event$Z[j],1))))
  return(a)
}
Event$D <- rep(0,nrow(Event))
for(j in 1:nrow(Event)){
  Event$D[j] <- Calculator(j)
}

rm(j)

Calculator <- function(j){
  a = (1-Event$X[j])*(ifelse(Event$Utype[j] %in% 1:4,0,
                             ifelse(Event$Utype[j] %in% 5:8,Event$D[j],
                                    ifelse(Event$Utype[j] %in% 9:12,1-Event$D[j],1)))) +
    (Event$X[j])*(ifelse(Event$Utype[j] %in% c(1,5,9,13),0,
                         ifelse(Event$Utype[j] %in% c(2,6,10,14),Event$D[j],
                                ifelse(Event$Utype[j] %in% c(3,7,11,15),1-Event$D[j],1))))
  return(a)
}
Event$Y <- rep(0,nrow(Event))
for(j in 1:nrow(Event)){
  Event$Y[j] <- Calculator(j)
}

rm(j)
# =================================================================================================== #
# Calculate p(D,X)
# =================================================================================================== #
p.10 = sum(Event$p[Event$Utype %in% 1:4 & Event$X == 0]) / sum(Event$p[Event$X == 0])
p.00 = 1 - (sum(Event$p[Event$Utype %in% 13:16 & Event$X == 0]) / sum(Event$p[Event$X == 0]))
p.11 = sum(Event$p[Event$Utype %in% c(1,5,9,13) & Event$X == 1]) / sum(Event$p[Event$X == 1])
p.01 = 1 - (sum(Event$p[Event$Utype %in% c(4,8,12,16) & Event$X == 1]) / sum(Event$p[Event$X == 1]))

ifelse(p.10 < p.11 & p.00 > p.01,1,0)
# =================================================================================================== #
# Calculate respective bounds
# =================================================================================================== #
Bound1 <- function(){
  Con1a = sum(Event$p[Event$Y == 0 & Event$X == 0 & Event$Z == 0]) /
    sum(Event$p[Event$X == 0 & Event$Z == 0])
  Con1b = sum(Event$p[Event$Y == 0 & Event$X == 0 & Event$Z == 1]) /
    sum(Event$p[Event$X == 0 & Event$Z == 1])
  Un1a = sum(Event$p[Event$Y == 0 & Event$Z == 0]) /
    sum(Event$p[Event$Z == 0])
  Un1b = sum(Event$p[Event$Y == 0 & Event$Z == 1]) /
    sum(Event$p[Event$Z == 1])
  
  a = ifelse(min(Con1a,Con1b) <= min(Un1a,Un1b),1,0)
  return(a)
}

Bound2 <- function(){
  Con1a = sum(Event$p[Event$Y == 0 & Event$D == 1 & Event$X == 1 & Event$Z == 0]) /
    sum(Event$p[Event$X == 1 & Event$Z == 0])
  Con1b = sum(Event$p[Event$Y == 0 & Event$D == 1 & Event$X == 1 & Event$Z == 1]) /
    sum(Event$p[Event$X == 1 & Event$Z == 1])
  Un1a = sum(Event$p[Event$Y == 0 & Event$D == 1 & Event$Z == 0]) /
    sum(Event$p[Event$Z == 0])
  Un1b = sum(Event$p[Event$Y == 0 & Event$D == 1 & Event$Z == 1]) /
    sum(Event$p[Event$Z == 1])
  
  a = ifelse(max(Con1a,Con1b) >= max(Un1a,Un1b),1,0)
  return(a)
}

Bound1()
Bound2()