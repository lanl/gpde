invlogit = function(x){
  exp(x)/(1+exp(x))
}

logit = function(p){
  if( any(p > 1 | p < 0) ) stop("Value of p must be in interval [0, 1].")
  log( p / (1-p) )
}
