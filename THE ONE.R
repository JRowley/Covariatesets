rm(list = ls())
m = 10000
p = 0.95
Source <- url("http://www.ucl.ac.uk/~zctpep9/Data%20archive/AE98Data.RData")
DF <- load(Source)															    # Load file from URL.
DF <- D																			        # Extract the data.
rm(list = c("Source", "D"))											    # Remove unnecessary objects from memory.
DF = data.frame(Y = DF$workedm, D = DF$morekids, X = DF$hispm, Z = DF$multi2nd)
n = nrow(DF)
library(iterpc)
library(Matrix)
library(mvtnorm)
library(plyr)
A = as.data.frame(table(DF))
A = ddply(A, .(X,Z), mutate, p = Freq / sum(Freq))

rm(list = c("DF"))

supp.D = length(unique(A$D))
supp.X = length(unique(A$X))
supp.Z = length(unique(A$Z))
supp.h = supp.D * supp.X
h = as.vector(outer(1:supp.D - 1, 1:supp.X - 1, paste, sep = ""))
h = as.vector(paste("h", h, sep = "."))

R = expand.grid(d = 1:supp.D - 1, x = 1:supp.X - 1)
S = expand.grid(xi = 1:supp.D - 1, eta = 1:supp.X - 1)

Q = expand.grid(eta = 1:supp.X - 1, z = 1:supp.Z - 1, d = 1:supp.D - 1, x = 1:supp.X - 1)

f.a <- function(q, r, s){
  eta = q[1]
  z = q[2]
  d = q[3]
  x = q[4]
  r.value = r$O[r$d %in% d & r$x %in% x]
  r.parameter = which(r$O == r.value)
  s.subset = s[s$eta %in% eta,]
  d.set = s.subset$xi[s.subset$O <= r.value]
  if(length(d.set) == 0){
    intermediate = c(w.event = rep(NA, nrow(A)), w.condition = rep(NA,nrow(A)))
  }
  if(length(d.set) != 0){
    w1 = (A$Y == 0)
    w2 = (A$X == eta)
    w3 = (A$Z == z)
    w4 = (A$D %in% d.set)
    w.event = w1 * w2 * w3 * w4
    w.condition = w2 * w3
    intermediate = c(w.event = w.event, w.condition = w.condition)
  }
  output = list(Info = c(h = r.parameter, sign = 1, O = r$O[r.parameter]), w = intermediate)
  return(output)
}
f.b <- function(q, r, s){
  eta = q[1]
  z = q[2]
  d = q[3]
  x = q[4]
  r.value = r$O[r$d %in% d & r$x %in% x]
  r.parameter = which(r$O == r.value)
  s.subset = s[s$eta %in% eta,]
  d.set = s.subset$xi[s.subset$O >= r.value]
  if(length(d.set) == 0){
    intermediate = c(w.event = rep(NA, nrow(A)), w.condition = rep(NA,nrow(A)))
  }
  if(length(d.set) != 0){
    w1 = (A$Y == 1)
    w2 = (A$X == eta)
    w3 = (A$Z == z)
    w4 = (A$D %in% d.set)
    w5 = 1 - (w1 * w4)
    w.event = w2 * w3 * w5
    w.condition = w2 * w3
    intermediate = c(w.event = w.event, w.condition = w.condition)
  }
  output = list(Info = c(h = r.parameter, sign = -1, O = r$O[r.parameter]), w = intermediate)
  return(output)
}
bounds <- function(O){
  r = cbind(R, O = O)
  s = cbind(S, O = O)
  intermediate.a = apply(Q, 1, f.a, r = r, s = s)
  intermediate.b = apply(Q, 1, f.b, r = r, s = s)
  Info.a = do.call("rbind", lapply(intermediate.a, "[[", 1))
  Info.b = do.call("rbind", lapply(intermediate.b, "[[", 1))
  w.a = do.call("rbind", lapply(intermediate.a, "[[", 2))
  w.b = do.call("rbind", lapply(intermediate.b, "[[", 2))
  Info = rbind(Info.a, Info.b)
  w = rbind(w.a, w.b)
  delete = which(is.na(apply(w, 1, max)) == 1)
  Info = data.frame(Info[-delete,])
  w = unname(w[-delete,])
  w.event = w[, 1:nrow(A)]
  w.condition = w[, (nrow(A) + 1):(2 * nrow(A))]
  Info$p = (w.event %*% A$p)
  output = list(Info = Info, w.event = w.event, w.condition = w.condition)
  return(output)
}

O = iterpc(n = 4, r = 4, ordered = TRUE)
O = getall(O)
colnames(O) = h
O = O[-c(7, 8, 13:24),]
Intermediate = apply(O, 1, bounds)

rep.TYPEA <- function(x){
  y = apply(x, 2, rep, nrow(x))
  return(y)
}
rep.TYPEB <- function(x){
  y = matrix(rep(c(x), each = nrow(x)), ncol = ncol(x))
  return(y)
}

builder <- function(o, g){
  j = which(apply(O, 1, function(x) all(x == o)) == 1)
  Input = Intermediate[[j]]
  Info = rbind(Input$Info, Input$Info)
  Info$Origin = c(rep(1, nrow(Input$Info)), rep(0, nrow(Input$Info)))
  delete = which(Info$p == 0)
  if(length(delete) > 0){
  Info = Info[-delete,]
  w.event = rbind(Input$w.event, Input$w.event)[-delete,]
  w.condition = rbind(Input$w.condition, Input$w.condition)[-delete,]
  }
  if(length(delete) == 0){
    Info = Info
    w.event = rbind(Input$w.event, Input$w.event)
    w.condition = rbind(Input$w.condition, Input$w.condition)
  }
  grid.subtraction = g[Info$h] * (1 - Info$Origin)
  Info$p = Info$p - grid.subtraction
  w.event.A = rep.TYPEA(w.event)
  w.event.B = rep.TYPEB(w.event)
  w.condition.A = rep.TYPEA(w.condition)
  w.condition.B = rep.TYPEB(w.condition)
  w.condition = w.condition.A * w.condition.B
  bound.1 = rep(Info$p, nrow(Info))
  bound.2 = rep(Info$p, each = nrow(Info))
  n.base = w.condition %*% A$Freq                                     # Denominator.
  n.base.1 = w.condition.A %*% A$Freq
  n.base.2 = w.condition.B %*% A$Freq
  n1 = (w.condition * w.event.A * w.event.B) %*% A$Freq               # Intersection.
  n2 = (w.condition * w.event.A * (1 - w.event.B)) %*% A$Freq         # A \ B.
  n3 = (w.condition * (1 - w.event.A) * w.event.B) %*% A$Freq         # B \ A.
  n4 = (w.condition * (1 - w.event.A) * (1 - w.event.B)) %*% A$Freq   # Complement.
  p1 = (1 - bound.1) * (1 - bound.2)
  p2 = (1 - bound.1) * (0 - bound.2)
  p3 = (0 - bound.1) * (1 - bound.2)
  p4 = (0 - bound.1) * (0 - bound.2)
  Variance = (n.base.1 ** -1) * (n1 * p1 + n2 * p2 + n3 * p3 + n4 * p4) * (n.base.2 ** -1)
  Variance = matrix(Variance, nrow = nrow(Info))
  output = list(Info = Info, Variance = Variance)
  return(output)
}
changer <- function(o, g){
  Input = builder(o, g)
  Info = Input$Info
  Variance = Input$Variance
  m1 = sapply(Info$sign, "<", Info$sign)
  m2 = sapply(Info$O, ">", Info$O)
  m3 = do.call(rbind, replicate(nrow(Info), Info$Origin, simplify = FALSE))
  m4 = do.call(cbind, replicate(nrow(Info), Info$Origin, simplify = FALSE))
  m = m1 * m2 * m3 * m4
  rm(m1, m2, m3, m4)
  m1 = bdiag(split(m, col(m)))
  m2 = diag(Info$Origin, nrow = nrow(Info), ncol = nrow(Info))
  m3 = do.call(rbind, replicate(nrow(Info), m2, simplify = FALSE))
  m4 = m3 - m1
  delete = which(rowSums(m1) == 0)
  m = m4[-delete,]
  rm(m1, m2, m3, m4, delete)
  m1 = m
  m2 = (1 - Info$Origin) * Info$sign
  m3 = diag(m2)
  delete = which(m2 == 0)
  m4 = m3[-delete,]
  m = rbind(m4, m1)
  rm(m1, m2, m3, m4)
  B = m %*% Info$p
  V = m %*% Variance %*% t(m)
  output = list(Estimate = B, Variance = V)
  return(output)
}
calculator <- function(B,V){
  b = as.matrix(B)
  v = as.matrix(V)
  s = sqrt(diag(v))
  normal = diag(s ** -1) %*% v %*% diag(s ** -1)
  Z = rmvnorm(m, mean = rep(0, nrow(v)), sigma = normal, method = "svd")
  K.aux = sort(apply(Z, 1, max))
  k.gamma = K.aux[floor(m * (1 - .1 / log(n)))]
  v.set = which((b - max(b - k.gamma * s) + 2 * k.gamma * s) >= 0)
  if (length(v.set == 1)){
    K = sort(Z[,v.set])
  }
  if (length(v.set > 1)){
    K = sort(apply(Z[,v.set], 1, max))
  }
  k = K[m * p]
  theta = max(b[v.set] - k * s[v.set])
  output = theta <= 0
  return(output)
}
engine <- function(o, g){
  Input = changer(o, g)
  B = Input$Estimate
  V = Input$Variance
  output = calculator(B, V)
  return(output)
}
final <- function(g){
  x = g
  output = apply(O, 1, engine, g = x)
  return(output)
}
system.time(final(c(0.01,0.01,0.01,0.01)))
