\name{kalman_ivp}
\alias{kalman_ivp}
\title{Probabilistic Initial Value Problem Solver Through Kalman Filter-Smoother}
\usage{
kalman_ivp((t,
            odefn,
            init_u,
            init_obs,
            n_draws,
            n_states   = 1,
            proc_var   = 0.2,
            proc_len   = dt*2,
            err        = 0.0,
            cov_type   = 'matern',
            obs_type   = 'mean',
            verbose    = FALSE))
}
\arguments{
  \item{t}{
    Vector of times for evaluation of the solver. Currently uses the first two times to get a fixed dt.
  }
  \item{odefn}{
    Function that takes two arguments, t and u. t is the current time and u is the vector state used to obtain the derivatives for the solver.
  }
  \item{init_u}{
    Initial state for the system to begin time integration.
  }
  \item{init_obs}{
    Boolean vector to indicate which states are observed at time zero.
  }
  \item{n_draws}{
    Currently no effect. Drawing multiple samples best done by wrapping this function in a form of apply statement.
  }
  \item{n_states}{
    Number of states in the ODE.
  }
  \item{proc_var}{
    Gaussian process variance for the probabilistic solver.
  }
  \item{proc_len}{
    Gaussian process correlaton length for the solver.
  }
  \item{err}{
    Process error term. Allows solutions to be imprecise at the cost of uncertainty in the solution.
  }
  \item{cov_type}{
    Covariance function for the Gaussian process of the probabilistic solver. matern and int_wiener are the only implemented options for the Kalman filter-based solver.
  }
  \item{obs_type}{
    Type of observation used for the solver. There are two options: "chkrebtii", which uses the approach of Chkrebtii 2013 by drawing a random value from the process to use for the observation or "mean", which uses the GP mean as the observation. The first fully incorporates solver uncertainty into the solution, while the second allows the mean solution to be consistent with rk4.
  }
  \item{verbose}{
    Boolean used to determine whether diagnostic information is printed.
  }
}
\description{
Evaluates a sample of the probabilistic solver using the Kalman filter formulation of a Gaussian process.
}
