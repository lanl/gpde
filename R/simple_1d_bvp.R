simple_1d_bvp = function(t,
                         ode,
                         bc,
                         proc_var   = 0.2,
                         proc_len   = 0.8,
                         cov_func   = "sq_exp",
                         verbose    = FALSE){

  if (cov_func == 'sq_exp') {
    cov_mat = proc_var * sqexp_cov_1d(t, proc_len)
  } else if (cov_func == 'abs_exp'){
    cov_mat = proc_var * absexp_cov_1d(t,proc_len)
  } else if (cov_func == 'matern32'){
    cov_mat = proc_var * matern_cov_1d(t, proc_len)
  } else{
    stop("Covariance Function Not Found")
  }

  y_ind = 1:n + n
  dyind = 1:n

  state_mean = matrix(c(rep(bc[1], n), rep(bc[2], n)), n, 2)
  deriv_mean = matrix(c(rep(    0, n), rep(bc[1], n)), n, 2)

  state_mean = cov_mat[y_ind, y_ind[n]] %*% solve(cov_mat[y_ind[n], y_ind[n]]) %*% c(bc[3], bc[4]) + state_mean
  deriv_mean = cov_mat[dyind, y_ind[n]] %*% solve(cov_mat[y_ind[n], y_ind[n]]) %*% c(bc[3], bc[4]) + deriv_mean
  cov_mat    = cov_mat - cov_mat[ ,y_ind[n]] %*% solve(cov_mat[y_ind[n], y_ind[n]]) %*% cov_mat[y_ind[n], ] + 1.e-8*diag(2*n)

  state_cov = cov_mat[y_ind, y_ind]
  deriv_cov = cov_mat[dyind, dyind]
  cross_cov = cov_mat[y_ind, dyind]

  gg       = deriv_cov[1,1]
  f_resid  = (odefn(t[1], bc) - c(0, bc[1]))/gg

  for (nn in 1:(n-1)){
    state_mean = state_mean + cross_cov[ ,nn] %o% f_resid
    deriv_mean = deriv_mean + deriv_cov[ ,nn] %o% f_resid

    state_cov = state_cov - (cross_cov[ ,nn]/gg) %*%  t(cross_cov[, nn])
    cross_cov = cross_cov - (cross_cov[ ,nn,drop=F]/gg) %*%  deriv_cov[nn, ,drop=F]
    deriv_cov = deriv_cov - (deriv_cov[ ,nn,drop=F]/gg) %*%  deriv_cov[nn, ,drop=F]

    ff = odefn(t[nn+1], state_mean[nn+1, ] + rnorm(2)*state_cov[nn+1, nn+1])
    gg = 2*deriv_cov[nn+1, nn+1]

    f_resid = (ff - deriv_mean[nn+1,])/gg
    if(verbose) cat(nn, " ")
  }
  if(verbose) cat("\n")
  return(sapply(1:n_draws,function(iii) t(chol(state_cov)) %*% matrix(rnorm(2*n),n,2) + state_mean))
  return()
}
