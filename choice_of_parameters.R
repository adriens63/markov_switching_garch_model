
library(MSGARCH)
library(xts)
simulate_MSGARCH <- function(n_periods , pi, A, mu, omega, alpha1, beta1) {
  
  # Initialize the state and volatility
  state <- rep(1,n_periods)
  state[1] <- sample(1:2, 1, prob = pi)
  sigma <- rep(0,n_periods)
  sigma[state] <- sqrt(omega[state[1]]/(1-alpha1[state[1]]-beta1[state[1]]))
  
  # Initialize the errors
  errors <- rnorm(n_periods, 0, 1)
  
  # Initialize the time series
  x <- rep(0,n_periods)
  x[1] <- mu[state[1]] + sigma[1] * errors[1]
  
  # Iterate through time
  for (i in 2:n_periods) {
    # Update the state
    state[i] <- sample(1:2, 1, prob = A[state[i-1],])
    
    # Update the volatility
    sigma[i] <- sqrt(omega[state[i]] + alpha1[state[i]] * (x[i-1]-mu[state[i-1]])^2 + beta1[state[i]] * sigma[i-1]^2)
    
    # Update the time series
    x[i] <- mu[state[i]] + sigma[i] * errors[i]
  }
  
  return(
    list(x = x,
         sigma = sigma,
         state = state,
         epsilon = errors)
  )
}


# Parameters
n_periods <- 50000
pi <- c(0.5, 0.5)
A <- matrix(c(0.98, 0.02, 0.04, 0.96), nrow = 2, byrow = TRUE) # transition matrix
mu <- c(0.06,-0.09)
omega <- c(0.3,2)
alpha1 <- c(0.35,0.10)
beta1 <- c(0.2, 0.6)


# Simulate the time series
simul_2states <- simulate_MSGARCH(n_periods = n_periods,
                                  pi = pi, 
                                  A = A,
                                  mu = mu,
                                  omega = omega,
                                  alpha1 = alpha1,
                                  beta1 = beta1)

library(ggplot2)
ggplot(data.frame(), aes(x= 1:1500, y= simul_2states$x[1:1500])) + geom_line() + xlab("Index") + ylab("Data Generating Process (1 to 1500 obs)")
plot(density(simul_2states$x),main = "Kernel density of simulated data")
acf(simul_2states$x^2,main= "ACF of squared data")

## Now we try to fit a 2 states and a 3-states ms-garch model to see whether there is a good fit
ms2.garch.spec <- CreateSpec(variance.spec = list(model = "sGARCH"),
                           distribution.spec = list(distribution = "norm"),
                           switch.spec = list(K = 2,do.mix = FALSE))

ms3.garch.spec <- CreateSpec(variance.spec = list(model = "sGARCH"),
                             distribution.spec = list(distribution = "norm"),
                             switch.spec = list(K = 3,do.mix = FALSE))
## fit without the bound constraints on the priors which seem to bias the optimization process
fit.2states <- FitMCMC(spec = ms2.garch.spec, data = simul_2states$x)
summary(fit.2states)

fit.3states <- FitMCMC(spec = ms3.garch.spec, data = simul_2states$x)
summary(fit.3states)



