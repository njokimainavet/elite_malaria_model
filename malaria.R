

library(tidyverse)
library(deSolve)

# -------- FUNCTION -------------
# 1. Define the SIR model
sir_si_itn <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    Nh  <-  Sh + Ih +Rh
    Nv  <-  Sv + Iv 
    dSh <-  -beta*Sh*(Nv/Nh)*(Iv/Nv)/4 + (alpha*Rh)
    dIh <-  beta*Sh*(Nv/Nh)*(Iv/Nv)/4 - (gammah*Ih)
    dRh <-  (gammah*Ih) - (alpha*Rh)
    dSv <-  r*Nv*(1-(Nv/K)) - beta*Sv*(Ih/Nh)/4 - (1+sin(t*2*pi/365))/5* d*Sv
    dIv <-  beta*Sv*(Ih/Nh)/4 - (1+sin(t*2*pi/365))/5*d*Iv
    
    return(list(c(dSh, dIh, dRh, dSv, dIv))) # , Nh, Nv
  })
}


# -------------------------
# 2. Define the parameters
# -------------------------
# Let's consider 5 years time period
dt <- (1:365)* 5 
R0 <- 2

# The infectious period is the average duration for which an individual remains infectious.
infectious_period <- 14
gammah <- 1/infectious_period

alpha <- 1/365

# The define beta, we will use the relation R0 = beta / gamma  because it is easier to interpret, that is, given an R0 and infectious period, we can calculate the corresponding beta value.
# beta is not directly interpretable because it is a composite of a number of factors.
beta = betah = betav <- R0 * gammah

K <- 5e5; d <- 1/28; r <- 0.08

# The parameters
params <- c(
  beta = beta,
  gamma = gammah,
  K = K, d = d, r = r
)

# ---------------------------------
# 3. Define the initial conditions
# ---------------------------------
Nh = 200000; Nv = 2*Nh
I0h =  1
I0v = 10000
# Initial conditions for S, I, R
inits <- c(Sh = Nh - I0h, Ih = I0h, Rh = 0,
           Sv = Nv - I0v, Iv = I0v)

results <- deSolve::lsoda(
  y = inits,
  times = dt,
  func = sir_si_itn,
  parms = params,
  hmax = 0.1        # <- make sure this is non-negative
)

results <- as.data.frame(results)

# Create data for ggplot2 by reshaping
results_long <- results |>
  pivot_longer(
    cols = c(2:6),
    names_to = "compartment",
    values_to = "value"
  )

sir_plot <- results_long %>% 
  ggplot(
  # data = 
  aes(
    x = time,
    y = value,
    color = compartment
  )
) +
  geom_line(linewidth = 1) +
  labs(
    title = "SIR model",
    x = "Time",
    y = "Number of individuals"
  )
plot(sir_plot)















