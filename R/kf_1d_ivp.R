kf_1d_ivp = function(t,
                     odefn,
                     init_u,
                     n_draws,
                     proc_var   = 0.2,
                     proc_len   = 0.8,
                     verbose    = FALSE,
                     step_plots = FALSE){
  n     = length(t)
  y_ind = 1:n + n
  dyind = 1:n


  cov_mat = proc_var*cov_1d(xx, proc_len)

  pred_u   = matrix(c(rep(init_u[1], n), rep(init_u[2], n)), n, 2)
  pred_du  = matrix(c(rep(        0, n), rep(init_u[1], n)), n, 2)

  state_cov = cov_mat[y_ind, y_ind]
  deriv_cov = cov_mat[dyind, dyind]
  cross_cov = cov_mat[y_ind, dyind]

  gg       = deriv_cov[1,1]
  f_resid  = (odefn(xx[1], init_u) - c(0, init_u[1]))/gg

  if(step_plots){
    plot_draws = sapply(1:10,function(iii) t(chol(state_cov)) %*% matrix(rnorm(2*n),n,2) + pred_u)
    par(mfrow = c(2,1))
    matplot(x    = xx,
            y    = plot_draws[1:n,],
            type = "l",
            lwd  = 2,
            ylab = "dy/dt")
    lines(x   = xx,
          y   = ode_du(xx),
          col = "purple",
          lwd = 3)
    matplot(x    = xx,
            y    = plot_draws[1:n + n,],
            type = "l",
            lwd  = 2,
            ylab = "y")
    lines(x   = xx,
          y   = ode_u(xx),
          col = "purple",
          lwd = 3)
  }

  for (nn in 1:(n-1)){
    pred_u  = pred_u  + cross_cov[ ,nn] %o% f_resid
    pred_du = pred_du + deriv_cov[ ,nn] %o% f_resid

    state_cov = state_cov - (cross_cov[ ,nn]/gg) %*%  t(cross_cov[, nn])
    cross_cov = cross_cov - (cross_cov[ ,nn,drop=F]/gg) %*%  deriv_cov[nn, ,drop=F]
    deriv_cov = deriv_cov - (deriv_cov[ ,nn,drop=F]/gg) %*%  deriv_cov[nn, ,drop=F]

    if(step_plots){
      plot_draws = sapply(1:10,function(iii) t(chol(state_cov)) %*% matrix(rnorm(2*n),n,2) + pred_u)
      par(mfrow  = c(2, 1))
      matplot(x    = xx,
              y    = plot_draws[1:n,],
              type = "l",
              lty  = 2,
              col  = "red",
              lwd  = 2,
              ylab = "dy/dt")
      lines(x   = xx,
            y   =  ode_du(xx),
            col = "purple",
            lwd = 3)
      matplot(x    = xx,
              y    = plot_draws[1:n + n,],
              type = "l",
              lty  = 2,
              col  = "red",
              lwd  = 2,
              ylab = "y")
      lines(x   = xx,
            y   = ode_u(xx),
            col = "purple",
            lwd = 3)
    }


    ff = odefn(xx[nn+1], pred_u[nn+1, ] + rnorm(2)*state_cov[nn+1, nn+1])
    gg = 2*deriv_cov[nn+1, nn+1]

    f_resid = (ff - pred_du[nn+1,])/gg
    if(verbose) cat(nn, " ")
  }
  if(verbose) cat("\n")
  return(sapply(1:n_draws,function(iii) t(chol(state_cov)) %*% matrix(rnorm(2*n),n,2) + pred_u))
}
