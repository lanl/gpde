# Re-tested 11/10/2022 and works on a basic level.

library(gpde)

n  = 100
xx = 0:(n - 1)/(n - 1)*10

init_u = c(0,-1)

odefn  = function(t,u) c(sin(2*t) - u[2], u[1])

plot_draws = simple_1d_ivp(t        = xx,
                           odefn    = odefn,
                           init_u   = init_u,
                           cov_func = 'sq_exp',
                           err      = 0.002,
                           proc_len = 0.5,
                           n_draws  = 200)


ode_u  = function(x) (-4*cos(x) + 2*sin(x) - sin(2*x) + cos(x))/3
ode_du = function(x) (4*sin(x) + 2*cos(x) - 2*cos(2*x) - sin(x))/3

par(mfrow  = c(2, 1))
par(mar    = c(4, 4,1,1))
matplot(x    = xx,
        y    = t(apply(plot_draws[1:n,],1,quantile,p=c(0.025,0.975))),
        type = "l", lty = 2,col="red",
        lwd  = 2,
        xlab = "X",
        ylab = "dy/dt")
lines(x   = xx,
      y   =  ode_du(xx),
      col = "purple",
      lwd = 3)
legend("bottomleft",
       legend = c("True Solution", "95% Credible Bounds"),
       lty    = 1:2,
       lwd    = 2,
       bty    = "n",
       col    = c("purple","red"))

matplot(x    = xx,
        y    = t(apply(plot_draws[1:n + n,],1,quantile,p=c(0.025,0.975))),
        type = "l", lty = 2,col="red",
        lwd  = 2,
        xlab = "X",
        ylab = "y")
lines(x   = xx,
      y   = ode_u(xx),
      col = "purple",
      lwd = 3)

