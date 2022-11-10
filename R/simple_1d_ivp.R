simple_1d_ivp = function(t,
                         odefn,
                         init_u,
                         n_draws,
                         proc_var   = 0.2,
                         proc_len   = 0.8,
                         cov_func   = "sq_exp",
                         err        = 0.0,
                         makeplots  = FALSE,
                         verbose    = FALSE){
  p     = length(init_u)
  n     = length(t)
  y_ind = 1:n + n
  dyind = 1:n


  if (cov_func == 'sq_exp') {
    cov_mat = proc_var * sqexp_cov_1d(t, proc_len)
  } else if (cov_func == 'abs_exp'){
    cov_mat = proc_var * absexp_cov_1d(t,proc_len)
  } else if (cov_func == 'matern32'){
    cov_mat = proc_var * matern_cov_1d(t, proc_len)
  } else{
    stop("Covariance Function Not Found")
  }

  # Set state means to initial boundary condition
  state_mean = lapply(1:n_draws,function(ii) sapply(init_u, function(x) rep(x,n)) )
  # Initialize the derivative mean as a zero mean GP
  deriv_mean = lapply(1:n_draws,function(ii) matrix(0, n, p))

  # Build the state and derivative covariance matrices and their cross covariance
  state_cov = cov_mat[y_ind, y_ind]
  deriv_cov = cov_mat[dyind, dyind]
  cross_cov = cov_mat[y_ind, dyind]

  gg       = deriv_cov[1,1]
  f_resid  = lapply(1:n_draws, function(ii) (odefn(t[1], init_u) - deriv_mean[[ii]][1,])/gg )

  for (nn in 1:(n-1)){
    state_mean = lapply(1:n_draws, function(ii) state_mean[[ii]] + cross_cov[ ,nn] %o% f_resid[[ii]] )
    deriv_mean = lapply(1:n_draws, function(ii) deriv_mean[[ii]] + deriv_cov[ ,nn] %o% f_resid[[ii]] )

    state_update = (cross_cov[ ,nn,drop=F]/gg) %*%  t(cross_cov[, nn,drop=F])
    cross_update = (cross_cov[ ,nn,drop=F]/gg) %*%  deriv_cov[nn, ,drop=F]
    deriv_update = (deriv_cov[ ,nn,drop=F]/gg) %*%  deriv_cov[nn, ,drop=F]

    state_cov = state_cov - state_update
    cross_cov = cross_cov - cross_update
    deriv_cov = deriv_cov - deriv_update

    if(makeplots){
      png(filename = paste0("frame_",nn,".png"),width = 1000,height = 700)
      plot_draws = sapply(state_mean,function(xx) t(chol(state_cov + 1.e-8*diag(n))) %*% matrix(rnorm(2*n),n,2) + xx)
      par(mfrow  = c(2, 1))
      matplot(x    = t,
              y    = plot_draws[1:n,],
              type = "l",
              lty  = 2,
              col  = rgb(1,0,0,0.3),
              lwd  = 2,
              ylab = "dy/dt",
              ylim = c(-4,4))
      lines(x   = t,
            y   =  odefn(t),
            col = "purple",
            lwd = 3)
      matplot(x    = t,
              y    = plot_draws[1:n + n,],
              type = "l",
              lty  = 2,
              col  = rgb(1,0,0,0.3),
              lwd  = 2,
              ylab = "y",
              ylim = c(-4,4))
      lines(x   = t,
            y   = ode_u(t),
            col = "purple",
            lwd = 3)
      dev.off()
    }

    ff = lapply(1:n_draws, function(ii) odefn(t[nn+1], state_mean[[ii]][nn+1, ] + rnorm(p)*sqrt(state_cov[nn+1, nn+1])) )
    gg = (1 + err)*deriv_cov[nn+1, nn+1]

    f_resid = lapply(1:n_draws, function(ii) (ff[[ii]] - deriv_mean[[ii]][nn+1,])/gg )
    if(verbose) cat(nn, " ")
  }
  if(verbose) cat("\n")
  for (jj in 1:1000){
    chol_mat = tryCatch({ chol(kronecker(diag(p),state_cov))}, error = function(cond) "Bad Value")
    if ( typeof(chol_mat) != "character") break
    #print("adding nugget")
    state_cov = state_cov + 1.e-6*diag(n)
  }
  # return(list(mean = state_mean, cov=state_cov))
  return( sapply(1:length(state_mean),
                 function(xx)
                   c(state_mean[[xx]]) + t(chol_mat)%*%rnorm(p*n) ) )
}
