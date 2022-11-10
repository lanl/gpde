# Re-tested 11/10/2022 and works at the basic level.
library(gpde)

# RK4 solver for the Duffing equation

rk4.duffing <- function(yini, par, h, max.t) {

  # Functions for derivatives

  du = function(t, u, v, par) {
    return(v)
  }

  dv = function(t, u, v, par) {
    d = par[1]
    gamma = par[2]
    omega = par[3]
    return(-d*v + u - u^3 + gamma*cos(omega*t))
  }

  # Create variables to store results

  t = seq(0, max.t, by = h)
  n = max.t/h + 1

  u = rep(0, n)
  v = rep(0, n)

  u[1] = yini[1]
  v[1] = yini[2]

  # Fourth order Runge-Kutta algorithm

  for(i in 1:(n-1)) {

    k11 = h*du(t[i], u[i], v[i], par)
    k12 = h*dv(t[i], u[i], v[i], par)

    k21 = h*du(t[i] + h/2, u[i] + k11/2, v[i] + k12/2, par)
    k22 = h*dv(t[i] + h/2, u[i] + k11/2, v[i] + k12/2, par)

    k31 = h*du(t[i] + h/2, u[i] + k21/2, v[i] + k22/2, par)
    k32 = h*dv(t[i] + h/2, u[i] + k21/2, v[i] + k22/2, par)

    k41 = h*du(t[i] + h, u[i] + k31, v[i] + k32, par)
    k42 = h*dv(t[i] + h, u[i] + k31, v[i] + k32, par)

    u[i + 1] = u[i] + 1/6*(k11 + 2*k21 + 2*k31 + k41)
    v[i + 1] = v[i] + 1/6*(k12 + 2*k22 + 2*k32 + k42)

  }

  return(list(t = t, u = u, v = v))

}

# Initial conditions, maximum time, and step size

h = 0.05
max.t = 1000

xx = seq(0, max.t, by = h)

n = length(xx)


init_u = c(0,0) # Initial conditions

d = 0.2
gamma = 0.3
omega = 1.0

odefn  = function(t,u) c(-d*u[1] + u[2] - u[2]^3 + gamma*cos(omega*t), u[1])

# Get fourth order Runge-Kutta solution

rk4.sol = rk4.duffing(init_u, c(d, gamma, omega),h/10,max.t)

# Get GP Solver solution
xx = seq(0, max.t, by = h/10)
kf_result = kalman_1d_ivp(t        = xx,
                          odefn    = odefn,
                          init_u   = init_u,
                          cov_type = 'int_wiener',
                          obs_type = 'mean',
                          err      = 0.001,
                          proc_len = 0.1,
                          proc_var = 5,
                          n_draws  = 2)

d_cov = apply(kf_result$cov, 3, function(x)x[2,2])
u_cov = apply(kf_result$cov, 3, function(x)x[1,1])
par(mfrow  = c(2, 1))
par(mar    = c(4, 4,1,1))
plot(x    = xx,
        y    = kf_result$mean[2,],
        type = "l", lty = 1,
        lwd  = 3,
        xlab = "t",
        ylab = "dy/dt", col = rgb(0,0,0,.5))
lines(xx, kf_result$mean[2,] + 2*sqrt(d_cov), type = 'l', col = 'blue')
lines(xx, kf_result$mean[2,] - 2*sqrt(d_cov), type = 'l', col = 'blue')
lines(xx, kf_result$mean[2,], type = 'l', col = 'green')
lines(rk4.sol$t, rk4.sol$v, type = 'l', col = 'red')
plot(x    = xx,
        y    = kf_result$mean[1,],
        type = "l", lty = 1,
        lwd  = 3,
        xlab = "t",
        ylab = "y", col = rgb(0,0,0,.5))
lines(xx, kf_result$mean[1,] + 2*sqrt(u_cov), type = 'l', col = 'blue')
lines(xx, kf_result$mean[1,] - 2*sqrt(u_cov), type = 'l', col = 'blue')
lines(rk4.sol$t, rk4.sol$u, type = 'l', col = 'red')
