Rate <- function(X){
  Conditional.Artstein = (2*X)*(2)*(factorial(2)*X)
  Unconditional.Artstein = (2*X)*(2)*(factorial(2*X))
  Report = Conditional.Artstein + Unconditional.Artstein
  return(Report)
}

Conditional.Rate <- function(X){
  Conditional.Artstein = (2*X)*(2)*(factorial(2)*X)
  return(Conditional.Artstein)
}