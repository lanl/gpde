kalman_1d_ivp = function(t,
                         odefn,
                         init_u,
                         n_draws,
                         proc_var   = NULL,
                         proc_len   = NULL,
                         err        = 0.0,
                         cov_type   = 'matern',
                         obs_type   = 'mean',
                         verbose    = FALSE){
  p     = length(init_u)
  n     = length(t)
  y_ind = 1:n + n
  dyind = 1:n
  dt = t[2] - t[1]

  if ( is.null(proc_var) ) proc_var = 0.2
  if ( is.null(proc_len) ) proc_len = dt*2

  # if (cov_func == 'sq_exp') {
  #   cov_mat = proc_var * sqexp_cov_1d(t, proc_len)
  # } else if (cov_func == 'abs_exp'){
  #   cov_mat = proc_var * absexp_cov_1d(t,proc_len)
  # } else if (cov_func == 'matern32'){
  #   cov_mat = proc_var * matern_cov_1d(t, proc_len)
  # } else{
  #   stop("Covariance Function Not Found")
  # }
  #
  # # Set state means to initial boundary condition
  # state_mean = lapply(1:n_draws,function(ii) sapply(init_u, function(x) rep(x,n)) )
  # # Initialize the derivative mean as a zero mean GP
  # deriv_mean = lapply(1:n_draws,function(ii) matrix(0, n, p))
  #
  # # Build the state and derivative covariance matrices and their cross covariance
  # state_cov = cov_mat[y_ind, y_ind]
  # deriv_cov = cov_mat[dyind, dyind]
  # cross_cov = cov_mat[y_ind, dyind]
  #
  # gg       = deriv_cov[1,1]
  # f_resid  = lapply(1:n_draws, function(ii) (odefn(t[1], init_u) - deriv_mean[[ii]][1,])/gg )

  lam   = sqrt(5)/proc_len
  if (cov_type == 'matern'){
    q          = proc_var/3/2*(2*lam)^(5)
    Q_mat      = matrix(0,3,3)
    Q_mat[3,3] = q*dt
    F_mat = diag(3) + matrix(c(0,1,0, 0,0,1, -lam^3,-3*lam^2,-3*lam), 3,3, byrow = TRUE) * dt
  } else if (cov_type == 'int_wiener'){
    F_mat = diag(3) + matrix(c(0,dt,dt^2,0,0,dt,0,0,0), 3, 3, byrow=TRUE)
    q          = proc_var
    Q_mat      = matrix(0,3,3)
    int_num    = 3
    for (ii in 1:3){
      for (jj in 1:3){
        F_mat[ii, jj] = ifelse(ii <= jj, dt^(jj - ii)/factorial(jj - ii), 0)
        Q_mat[ii, jj] = dt^(2*int_num + 1 - ii - jj)/(2*int_num + 1 - ii - jj)/factorial(int_num - ii)/factorial(int_num - jj)
      }
    }
    Q_mat = Q_mat*q
  } else{
    stop("cov_type not found. 'matern' and 'int_wiener' are only available options")
  }
  L_mat = matrix(c(0,0,1), 3, 1)
  H_mat = matrix(1,3,1)
  S_mat = diag(3)*err


  state_mean = matrix(0,3,n)
  state_cov  = array(0,dim=c(3,3,n))

  state_mean[ , 1]  = c(init_u[2], odefn(t[1], init_u)[2:1])
  state_cov[ , , 1] = diag(3)*1.e-8

  for (nn in 1:(n-1)){
    m_tmp = F_mat %*% state_mean[ , nn, drop = F]
    c_tmp = F_mat %*% state_cov[ , , nn] %*% t(F_mat) + Q_mat

    if (obs_type == 'chkrebtii'){
      y_draw = c( m_tmp[1:2,,drop=F] + t(chol(c_tmp[-3,-3])) %*% rnorm(2) )
      y_obs  = matrix(c(y_draw[1], odefn(t[nn+1], y_draw[2:1])[2:1]), 3, 1)
    } else if (obs_type == 'mean'){
      y_obs  = matrix(c(m_tmp[1,,drop=F], odefn(t[nn+1], m_tmp[2:1,,drop=F])[2:1]), 3, 1)
    } else {
      stop('Observation Type Not Found. Set "chkrebtii" or "mean".')
    }

    f_mat  = c_tmp + S_mat

    state_mean[ , nn+1]  = m_tmp + c_tmp %*% ( solve(f_mat) %*% (y_obs - m_tmp) )
    state_cov[ , , nn+1] = c_tmp - c_tmp %*% ( solve(f_mat) %*% c_tmp )
  }
  for (nn in (n-1):1){
    m_tmp = F_mat %*% state_mean[ , nn, drop = F]
    c_tmp = F_mat %*% state_cov[ , , nn] %*% t(F_mat) + Q_mat
    G_mat = state_cov[ , , nn] %*% t(F_mat) %*% solve(c_tmp)
    state_mean[ , nn]  = state_mean[ , nn] + G_mat %*% (state_mean[ , nn+1] - m_tmp)
    state_cov[ , , nn] = state_cov[ , , nn] - G_mat %*% (state_cov[ , , nn +1 ] - c_tmp) %*% t(G_mat)
  }
  if(verbose) cat("\n")
  # for (jj in 1:1000){
  #   chol_mat = tryCatch({ chol(kronecker(diag(p),state_cov))}, error = function(cond) "Bad Value")
  #   if ( typeof(chol_mat) != "character") break
  #   print("adding nugget")
  #   state_cov = state_cov + 1.e-6*diag(n)
  # }
  return(list(mean = state_mean, cov=state_cov))
  return( sapply(1:length(state_mean),
                 function(xx)
                   c(state_mean[[xx]]) + t(chol_mat)%*%rnorm(p*n) ) )
}
