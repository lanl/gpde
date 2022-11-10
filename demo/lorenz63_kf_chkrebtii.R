# Re-tested 11/10/2022 and works at basic level

library(gpde)

n        = 2001
xx       = 0:(n - 1)/(n - 1)*20
alpha    = 1/n
n_states = 3

odefn  = function(t,u) {
  sigma = 10
  rho   = 28
  beta  = 8/3
  ode_out = list()
  ode_out[[1]] = matrix(c(NA,sigma*(u[[2]][1] - u[[1]][1]),NA),3,1)
  ode_out[[2]] = matrix(c(NA,u[[1]][1]*(rho - u[[3]][1]) - u[[2]][1],NA),3,1)
  ode_out[[3]] = matrix(c(NA,u[[1]][1]*u[[2]][1] - beta*u[[3]][1],NA),3,1)
  return(ode_out)
}

init_mean = c(-12, -5, 38)
init_u    = odefn(0, lapply(init_mean,function(x) x))

for (ii in 1:n_states) {
  init_u[[ii]][1] = init_mean[ii]
  init_u[[ii]][3] = 0
}

get_solution = function(iii){
  tmp_sol = kalman_ivp(t        = xx,
                       odefn    = odefn,
                       init_u   = init_u,
                       init_obs = c(T,T,F),
                       cov_type = 'matern',
                       obs_type = 'chkrebtii',
                       n_states = 3,
                       err      = 0.00,
                       proc_var = alpha,
                       proc_len = 4*(xx[2] - xx[1]),
                       n_draws  = 50)
  return( cbind(tmp_sol$mean[[1]][1,],
                tmp_sol$mean[[2]][1,],
                tmp_sol$mean[[3]][1,]) )
}

plot_draws = sapply(1:25, get_solution)

par(mfrow  = c(3, 1))
par(mar    = c(4, 4,1,1))
matplot(x    = xx,
        y    = plot_draws[1:n,],
        type = "l", lty = 1,col=rgb(0,0,0,0.2),
        lwd  = 2,
        xlab = "t",
        ylab = "X")
matplot(x    = xx,
        y    = plot_draws[1:n + n,],
        type = "l", lty = 1,col=rgb(0,0,0,0.2),
        lwd  = 2,
        xlab = "t",
        ylab = "Y")
matplot(x    = xx,
        y    = plot_draws[1:n + 2*n,],
        type = "l", lty = 1,col=rgb(0,0,0,0.2),
        lwd  = 2,
        xlab = "t",
        ylab = "Z")

par(mfrow=c(1,1))
plot(x    = plot_draws[1:n, 1],
     y    = plot_draws[1:n + n, 1],
     type = 'l',
     col  = rgb(0,0,0,0.1),
     lwd  = 2)
for (jj in 2:25) lines(x    = plot_draws[1:n, jj],
                        y    = plot_draws[1:n + n, jj],
                        type = 'l',
                        col  = rgb(0,0,0,0.1),
                        lwd  = 2)

par(mfrow=c(1,1))
plot(x    = plot_draws[1:n, 1],
     y    = plot_draws[1:n + 2*n, 1],
     type = 'l',
     col  = rgb(0,0,0,0.1),
     lwd  = 2)
for (jj in 2:25) lines(x    = plot_draws[1:n, jj],
                       y    = plot_draws[1:n + 2*n, jj],
                       type = 'l',
                       col  = rgb(0,0,0,0.1),
                       lwd  = 2)

par(mfrow=c(1,1))
plot(x    = plot_draws[1:1000, 1],
     y    = plot_draws[1:1000 + n, 1],
     type = 'l',
     col  = rgb(0,0,0,0.1),
     lwd  = 2)
for (jj in 2:25) lines(x    = plot_draws[1:1000, jj],
                       y    = plot_draws[1:1000 + n, jj],
                       type = 'l',
                       col  = rgb(0,0,0,0.1),
                       lwd  = 2)

par(mfrow=c(1,1))
par(mar    = c(4, 4,4,3))
plot(x    = plot_draws[500:1200, 1],
     y    = plot_draws[500:1200 + 2*n, 1],
     type = 'l',
     col  = rgb(0,0,0,0.1),
     lwd  = 2,
     xlab = "X",
     ylab = "Z",
     main = "Lorenz Attractor")
for (jj in 2:25) lines(x    = plot_draws[500:1200, jj],
                       y    = plot_draws[500:1200 + 2*n, jj],
                       type = 'l',
                       col  = rgb(0,0,0,0.1),
                       lwd  = 2)
