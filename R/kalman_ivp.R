kalman_ivp = function(t,
                      odefn,
                      init_u,
                      init_obs,
                      n_draws,
                      n_states   = 1,
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


  state_mean = lapply(1:n_states,function(x) matrix(0,3,n) )
  state_cov  = lapply(1:n_states,function(x) array(0,dim=c(3,3,n)) )

  for (ii in 1:n_states) {
    state_mean[[ii]][ , 1]  = init_u[[ii]]
    state_cov[[ii]][ , , 1] = diag(3)*1.e-8 + diag((!init_obs)*1000)
  }

  for (nn in 1:(n-1)){
    m_tmp = lapply(state_mean, function(XX) F_mat %*% XX[ , nn, drop = F] )
    c_tmp = lapply(state_cov,  function(XX) F_mat %*% XX[ , , nn] %*% t(F_mat) + Q_mat )

    if (obs_type == 'chkrebtii'){
      y_draw = lapply(1:n_states, function(iii) c( m_tmp[[iii]] + t(chol(c_tmp[[iii]])) %*% rnorm(nrow(c_tmp[[iii]])) ) )
      y_obs  = odefn(t[nn+1], y_draw)
    } else if (obs_type == 'mean'){
      y_obs  = odefn(t[nn+1], m_tmp)
    } else {
      stop('Observation Type Not Found. Set "chkrebtii" or "mean".')
    }

    f_mat  = lapply(c_tmp, function(XX) XX + S_mat)

    for (ii in 1:n_states) {
      obs_inds = !is.na(y_obs[[ii]])
      state_mean[[ii]][ , nn+1]  = m_tmp[[ii]] + c_tmp[[ii]][,obs_inds] %*% ( solve(f_mat[[ii]][obs_inds,obs_inds]) %*% (y_obs[[ii]][obs_inds] - m_tmp[[ii]][obs_inds]) )
      state_cov[[ii]][ , , nn+1] = c_tmp[[ii]] - c_tmp[[ii]][,obs_inds] %*% ( solve(f_mat[[ii]][obs_inds,obs_inds]) %*% c_tmp[[ii]][obs_inds,] )
    }
  }
  for (nn in (n-1):1){
    m_tmp = lapply(state_mean, function(XX) F_mat %*% XX[ , nn, drop = F] )
    c_tmp = lapply(state_cov,  function(XX) F_mat %*% XX[ , , nn] %*% t(F_mat) + Q_mat)
    for (ii in 1:n_states){
      G_mat = state_cov[[ii]][ , , nn] %*% t(F_mat) %*% solve(c_tmp[[ii]])
      state_mean[[ii]][ , nn]  = state_mean[[ii]][ , nn]  + G_mat %*% (state_mean[[ii]][ , nn+1] - m_tmp[[ii]])
      state_cov[[ii]][ , , nn] = state_cov[[ii]][ , , nn] - G_mat %*% (state_cov[[ii]][ , , nn +1 ] - c_tmp[[ii]]) %*% t(G_mat)
    }
  }
  if(verbose) cat("\n")

  return(list(mean = state_mean, cov=state_cov))
  return( sapply(1:length(state_mean),
                 function(xx)
                   c(state_mean[[xx]]) + t(chol_mat)%*%rnorm(p*n) ) )
}
