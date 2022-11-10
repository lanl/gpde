erf = function(x) 2 * pnorm(x * sqrt(2)) - 1

sqexp_cov_1d = function(x,beta){
  n       = length(x)
  cov_mat = diag(2*n)

  for (ii in 2:n){
    for (jj in 1:(ii-1)){
      cov_mat[ii, jj] = sqrt(pi) * beta * exp( -(x[ii] - x[jj])^2/(4*beta^2) )
      cov_mat[jj, ii] = cov_mat[ii, jj]
    }
  }

  for (ii in 1:n){
    for (jj in 1:n){
      cov_mat[ii, jj + n] = pi * beta^2 * ( erf( (x[jj] - x[ii])/(2*beta) ) + erf( (x[ii] - 0)/(2*beta) ) )
    }
  }

  cov_mat[1:n + n, 1:n] = t( cov_mat[1:n, 1:n + n] )

  for (ii in 1:n){
    for (jj in 1:ii){
      term1 = pi*beta^2 * ( x[ii]*erf(x[ii]/beta/2) - (x[jj] - x[ii])*erf( (x[jj] - x[ii])/beta/2 ) + x[jj]*erf(x[jj]/beta/2))
      term2 = 2*sqrt(pi)*beta^3 * ( exp(-(x[ii]/beta/2)^2) - exp(-(x[jj] - x[ii])^2/beta^2/4) + exp(-(x[jj]/beta/2)^2) - 1)
      cov_mat[ii + n, jj + n] = term1 + term2
      cov_mat[jj + n, ii + n] = cov_mat[ii + n, jj + n]
    }
  }

  diag(cov_mat) = diag(cov_mat)

  cov_mat
}

# old_sqexp_cov_1d = function(x,beta){
#   n       = length(x)
#   cov_mat = diag(2*n)
#
#   for (ii in 2:n){
#     for (jj in 1:(ii-1)){
#       cov_mat[ii, jj] = sqrt(2*pi) * beta * dnorm(x[ii], x[jj], beta)
#       cov_mat[jj, ii] = cov_mat[ii, jj]
#     }
#   }
#
#   for (ii in 1:n){
#     for (jj in 1:n){
#       cov_mat[ii, jj + n] = sqrt(2*pi) * beta * ( pnorm(x[jj], x[ii], beta) - pnorm(0, x[ii], beta) )
#     }
#   }
#
#   cov_mat[1:n + n,1:n] = t( cov_mat[1:n, 1:n + n] )
#
#   for (ii in 1:n){
#     for (jj in 1:ii){
#       cov_mat[ii+ n, jj + n] = sqrt(2*pi) * beta * integrate( f     = function(z) pnorm(x[jj], z, beta) - pnorm(0, z, beta),
#                                                               lower = 0,
#                                                               upper = x[ii])$value
#       cov_mat[jj + n, ii + n] = cov_mat[ii + n, jj + n]
#     }
#   }
#
#   diag(cov_mat) = diag(cov_mat)
#
#   cov_mat
# }

mat32_cov_func = function(x2,x1,beta){
  d = abs(x1 - x2)
  (1 + sqrt(3)*d/beta)*exp(-sqrt(3)*d/beta)
}

matern_cov_1d = function(x,beta){
  n       = length(x)
  cov_mat = diag(2*n)

  for (ii in 2:n){
    for (jj in 1:(ii-1)){
      d = abs(x[ii] - x[jj])
      cov_mat[ii, jj] = (1 + sqrt(3)*d/beta)*exp(-sqrt(3)*d/beta)
      cov_mat[jj, ii] = cov_mat[ii, jj]
    }
  }

  for (ii in 1:n){
    for (jj in 1:n){
      if (x[ii] >= x[jj]){
        term1 = exp(-sqrt(3)*(x[ii] - x[jj])/beta) * (2*sqrt(3)*beta + 3*(x[ii] - x[jj]))
        term2 = exp(-sqrt(3)*x[ii]/beta)*(2*sqrt(3)*beta + 3*x[ii])
        cov_mat[ii, jj + n] = 1/3*(term1 - term2)
      } else {
        term1 = 2*beta/sqrt(3) - 1/3*exp(-sqrt(3)*x[ii]/beta)*(2*sqrt(3)*beta + 3*x[ii])
        term2 = 2*beta/sqrt(3) - 1/3*exp(-sqrt(3)*(x[jj] - x[ii])/beta)*(2*sqrt(3)*beta + 3*(x[jj] - x[ii]))
        cov_mat[ii, jj + n] = term1 + term2
      }
    }
  }

  cov_mat[1:n + n,1:n] = t( cov_mat[1:n, 1:n + n] )

  for (ii in 1:n){
    for (jj in 1:ii){
      if (x[ii] >= x[jj]){
        term1 = 1/3*beta*(exp(-sqrt(3)*x[jj]/beta)*(3*beta + sqrt(3)*x[jj]) - 3*beta + 2*sqrt(3)*x[jj])
        term2 = 2*sqrt(3)/3*beta*x[jj] + sqrt(3)/3*beta*(x[jj]-x[ii])*exp(sqrt(3)/beta*(x[jj]-x[ii])) + x[ii]*sqrt(3)/3*beta*exp(-sqrt(3)/beta*x[ii]) - beta^2*exp(sqrt(3)/beta*(x[jj]-x[ii])) + beta^2*exp(-sqrt(3)/beta*x[ii])
        cov_mat[ii + n, jj + n] = term1 + term2
        cov_mat[jj + n, ii + n] = cov_mat[ii + n, jj + n]
      } else {
        term1 = 1/3*beta*(exp(-sqrt(3)*x[ii]/beta)*(3*beta + sqrt(3)*x[ii]) - 3*beta + 2*sqrt(3)*x[ii])
        term2 = 2*sqrt(3)/3*beta*x[ii] + sqrt(3)/3*beta*(x[ii]-x[jj])*exp(sqrt(3)/beta*(x[ii]-x[jj])) + x[jj]*sqrt(3)/3*beta*exp(-sqrt(3)/beta*x[jj]) - beta^2*exp(sqrt(3)/beta*(x[ii]-x[jj])) + beta^2*exp(-sqrt(3)/beta*x[jj])
        cov_mat[ii + n, jj + n] = term1 + term2
        cov_mat[jj + n, ii + n] = cov_mat[ii + n, jj + n]

      }
    }
  }

  diag(cov_mat) = diag(cov_mat)

  cov_mat
}

absexp_cov_1d = function(x,beta){
  n       = length(x)
  cov_mat = diag(2*n)

  for (ii in 2:n){
    for (jj in 1:(ii-1)){
      d = abs(x[ii] - x[jj])
      cov_mat[ii, jj] = exp(-d/beta)
      cov_mat[jj, ii] = cov_mat[ii, jj]
    }
  }

  for (ii in 1:n){
    for (jj in 1:n){
      if (x[jj] <= x[ii]){
        cov_mat[ii, jj + n] = beta*exp(1/beta*(x[jj]-x[ii])) - beta*exp(-x[ii]/beta)
      } else {
        cov_mat[ii, jj + n] = beta*( 1 - exp(-x[ii]/beta) ) + beta*(1 - exp((x[ii] - x[jj])/beta))
      }
    }
  }

  cov_mat[1:n + n,1:n] = t( cov_mat[1:n, 1:n + n] )

  for (ii in 1:n){
    for (jj in 1:ii){
      if (x[ii] >= x[jj]){
        cov_mat[ii+ n, jj + n]  = beta*(beta*(exp(-x[jj]/beta) - 1) + x[jj]) + x[jj]*beta - beta^2*exp(1/beta*(x[jj]-x[ii])) + beta^2*exp(-x[ii]/beta)
        cov_mat[jj + n, ii + n] = cov_mat[ii + n, jj + n]
      } else {
        cov_mat[ii+ n, jj + n]  = beta*(beta*(exp(-x[ii]/beta) - 1) + x[ii]) + x[ii]*beta - beta^2*exp(1/beta*(x[ii]-x[jj])) + beta^2*exp(-x[jj]/beta)
        cov_mat[jj + n, ii + n] = cov_mat[ii + n, jj + n]

      }
    }
  }

  diag(cov_mat) = diag(cov_mat)

  cov_mat
}

