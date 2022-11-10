# Re-tested 11/10/2022 and works at a basic level.
library(gpde)

n  = 501
xx = 0:(n - 1)/(n - 1)*10
alpha = 1/n

init_u = c(-12, -5, 38)

odefn  = function(t,u) {
  sigma = 10
  rho   = 28
  beta  = 8/3
  return(
    c(sigma*(u[2] - u[1]),
    u[1]*(rho - u[3]) - u[2],
    u[1]*u[2] - beta*u[3])
  )
}

plot_draws = simple_1d_ivp(t        = xx,
                           odefn    = odefn,
                           init_u   = init_u,
                           cov_func = 'sq_exp',
                           err      = 0.00,
                           proc_var = alpha,
                           proc_len = 4*(xx[2] - xx[1]),
                           n_draws  = 10)

ode_u  = function(x) (-4*cos(x) + 2*sin(x) - sin(2*x) + cos(x))/3
ode_du = function(x) (4*sin(x) + 2*cos(x) - 2*cos(2*x) - sin(x))/3

par(mfrow  = c(3, 1))
par(mar    = c(4, 4,1,1))
matplot(x    = xx,
        y    = plot_draws[1:n,],
        type = "l", lty = 1,col="red",
        lwd  = 1,
        xlab = "t",
        ylab = "X")
matplot(x    = xx,
        y    = plot_draws[1:n + n,],
        type = "l", lty = 1,col="red",
        lwd  = 1,
        xlab = "t",
        ylab = "Y")
matplot(x    = xx,
        y    = plot_draws[1:n + 2*n,],
        type = "l", lty = 1,col="red",
        lwd  = 1,
        xlab = "t",
        ylab = "Z")

