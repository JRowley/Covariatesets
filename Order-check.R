rm(list = ls())
Source <- url("http://www.ucl.ac.uk/~zctpep9/Data%20archive/AE98Data.RData")
DF <- load(Source)															    # Load file from URL.
DF <- D																			        # Extract the data.
rm(list = c("Source","D"))											    # Remove unnecessary objects from memory.
DF = data.frame(Y = DF$workedm, D = DF$morekids, X = DF$hispm, Z = DF$multi2nd)
n = nrow(DF)
library(plyr)
library(Matrix)
library(mvtnorm)
A = as.data.frame(table(DF))
A = ddply(A, .(X,Z), mutate, p = Freq / sum(Freq))

rm(list = c("DF"))

supp.D = length(unique(A$D))
supp.X = length(unique(A$X))
supp.Z = length(unique(A$Z))
supp.h = supp.D * supp.X
h = as.vector(outer(1:supp.D - 1, 1:supp.X - 1, paste, sep = ""))
h = as.vector(paste("h", h, sep = "."))

bounds <- function(O){
  r = expand.grid(d = 1:supp.D - 1, x = 1:supp.X - 1)
  r$O = O
  s = r
  colnames(s) <- c("a","eta","O")
  f.a <- function(y){
    eta = y[1]
    z = y[2]
    d = y[3]
    x = y[4]
    r.sub = r$O[r$d == d & r$x == x]
    r.int = which(r$O == r.sub)
    s.sub = s[s$eta == eta,]
    s.int = s.sub$a[s.sub$O <= r.sub]
    if(length(s.int) == 0){
      intermediate.1 = list(h = r.int, p = NA, n = NA, n.c = NA)
      intermediate.2 = c(w = rep(NA, nrow(A)), w.prime = rep(NA,nrow(A)))
      output = list(intermediate.1, intermediate.2)
    }
    if(length(s.int) != 0){
      w1 = (A$Y == 0)
      w2 = (A$X == eta)
      w3 = (A$Z == z)
      w4 = (A$D %in% s.int)
      w = w1 * w2 * w3 * w4
      w.prime = w2 * w3
      intermediate.1 = list(h = r.int, p = sum(A$p * w), n = sum(A$Freq * w), n.c = sum(A$Freq * w.prime) - sum(A$Freq * w))
      intermediate.2 = c(w = w, w.prime = w.prime)
      output = list(intermediate.1, intermediate.2)
    }
    return(output)
  }
  f.b <- function(y){
    eta = y[1]
    z = y[2]
    d = y[3]
    x = y[4]
    r.sub = r$O[r$d == d & r$x == x]
    r.int = which(r$O == r.sub)
    s.sub = s[s$eta == eta,]
    s.int = s.sub$a[s.sub$O >= r.sub]
    if(length(s.int) == 0){
      intermediate.1 = list(h = r.int, p = NA, n = NA, n.c = NA)
      intermediate.2 = c(w = rep(NA, nrow(A)), w.prime = rep(NA,nrow(A)))
      output = list(intermediate.1, intermediate.2)
    }
    if(length(s.int) != 0){
      w1 = (A$Y == 1)
      w2 = (A$X == eta)
      w3 = (A$Z == z)
      w4 = (A$D %in% s.int)
      w5 = 1 - w1 * w4
      w = w2 * w3 * w5
      w.prime = w2 * w3
      intermediate.1 = list(h = r.int, p = sum(A$p * w), n = sum(A$Freq * w), n.c = sum(A$Freq * w.prime) - sum(A$Freq * w))
      intermediate.2 = c(w = w, w.prime = w.prime)
      output = list(intermediate.1, intermediate.2)
    }
    return(output)
  }
  q = expand.grid(eta = 1:supp.X - 1, z = 1:supp.Z - 1, d = 1:supp.D - 1, x = 1:supp.X - 1)
  int.a = apply(q, 1, f.a)
  int.b = apply(q, 1, f.b)
  inter.a = do.call("rbind", lapply(int.a, "[[", 1))
  inter.b = do.call("rbind", lapply(int.b, "[[", 1))
  inter.a = cbind(inter.a, sign = 1)
  inter.b = cbind(inter.b, sign = -1)
  int.prob = rbind(inter.a, inter.b)
  int.prob = int.prob[which(1 - apply(is.na(int.prob), 1, max) == 1),]
  int.prob = data.frame(int.prob)
  # The inequalities...
  i.a = do.call("rbind", lapply(int.a, "[[", 2))
  i.b = do.call("rbind", lapply(int.b, "[[", 2))
  i = rbind(i.a, i.b)
  i = i[which(1 - apply(is.na(i), 1, max) == 1),]
  i.w = i[,1:nrow(A)]
  i.prime = i[,(nrow(A) + 1):32]
  # Want to define four events: intersection, one occurs other doesn't, other occurs, neither occur.
  exp.w.1 = matrix(rep(c(i.w), each = nrow(i.w)), ncol = ncol(i.w))
  exp.w.2 = apply(i.w, 2, rep, nrow(i.w))
  exp.prime.1 = matrix(rep(c(i.prime), each = nrow(i.prime)), ncol = ncol(i.prime))
  exp.prime.2 = apply(i.prime, 2, rep, nrow(i.prime))
  E.base = exp.prime.1 * exp.prime.2
  E.inter.w = (E.base * exp.w.1 * exp.w.2) %*% A$Freq
  E.m2.w = (E.base * (exp.w.1 - exp.w.1 * exp.w.2)) %*% A$Freq
  E.m1.w = (E.base * (exp.w.2 - exp.w.1 * exp.w.2)) %*% A$Freq
  E.none.w = (E.base * (1 - exp.w.1) * (1 - exp.w.2)) %*% A$Freq
  p.1 = rep(as.numeric(int.prob$p), each = nrow(int.prob))
  p.2 = rep(as.numeric(int.prob$p), nrow(int.prob))
  variance = ((exp.prime.1 %*% A$Freq) ** -1) * (E.inter.w * (1 - p.1) * (1 - p.2) + E.m2.w * (1 - p.1) * (0 - p.2) + E.m1.w * (0 - p.1) * (1 - p.2) + E.none.w * p.1 * p.2) * ((exp.prime.2 %*% A$Freq) ** -1)
  sgn.1 = rep(as.numeric(int.prob$sign), each = nrow(int.prob))
  sgn.2 = rep(as.numeric(int.prob$sign), nrow(int.prob))
  sgn = sgn.1 * sgn.2
  variance = variance * sgn
  variance = matrix(variance, nrow = nrow(int.prob), byrow = FALSE)
  delete = which(apply(variance, 1, sum) == 0)
  int.prob = int.prob[-delete,]
  variance = variance[-delete, -delete]
  output = list(estimates = int.prob, variance = variance)
  return(output)
}
builder <- function(O){
  x = bounds(O)$estimates
  o = O[as.numeric(x$h)]
  x = cbind(x, O = o)
  y = bounds(O)$variance
  w1 = sapply(x$sign, "<", x$sign)
  w2 = sapply(x$O, ">=", x$O)
  w.prime = bdiag(split(w1 * w2, col(w1 * w2)))
  w3 = diag(1, nrow = nrow(x), ncol = nrow(x))
  w3 = do.call(rbind, replicate(nrow(x), w3, simplify = FALSE))
  w = w3 - w.prime
  delete = which(apply(w.prime, 1, sum) == 0)
  w = w[-delete,]
  B = w %*% as.numeric(x$p)
  V = w %*% y %*% t(w)
  output = list(estimates = B, variance = V)
  return(output)
}

O = iterpc(n = 4, r = 4, ordered = TRUE)
O = getall(O)
Results = apply(O, 1, builder)

m = 10000
p = 0.95

calculator <- function(B,V){
  b = as.matrix(B)
  v = as.matrix(V)
  s = sqrt(diag(V))
  normal = diag(s ** -1) %*% v %*% diag(s ** -1)
  Z = rmvnorm(m, mean = rep(0, nrow(v)), sigma = normal, method = "svd")
  K.aux = sort(apply(Z, 1, max))
  k.gamma = K.aux[floor(m * (1 - .1 / log(n)))]
  v.set = which((B - max(B - k.gamma * s) + 2 * k.gamma * s) >= 0)
  K = sort(apply(Z[,v.set], 1, max))
  k = K[m * p]
  theta = max(B[v.set] - k * s[v.set])
  output = theta <= 0
  return(output)
}
f <- function(Results){
  output = vector(length = length(Results))
  for(j in 1:length(Results)){
    B = Results[[j]][[1]]
    V = Results[[j]][[2]]
    output[j] = calculator(B,V)
  }
  return(output)
}
f(Results)