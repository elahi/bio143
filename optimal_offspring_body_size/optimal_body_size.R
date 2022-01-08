################################################################################
##' @title Optimal body size
##' @author Robin Elahi
##' @contact elahi@stanford.edu
##' @date 2022-01-07
##' @log 
################################################################################

##' Exercise from the Stanford course Quantitative Methods for Marine Ecology and Conservation (BIOHOPK143), 
##' taught by Giulio DeLeo, Maurice Goodman, and Richard Grewelle. 

##### PACKAGES, DATA #####
library(tidyverse)
library(here)

theme_set(theme_bw(base_size = 10) + 
            theme(panel.grid = element_blank(), 
                  strip.background = element_blank())
)

##### GRAPHICAL SOLUTION #####

# Set model parameters
delta <- 30 # total biomass available for reproduction (per adult)
gamma <- 1 # minimum offspring size necessary for survival
alpha <- 0.1 # parameter
beta  <- 0.25 # parameter

## Example: For body size x, we can compute lambda (finite growth rate)

# First, compute the fraction of offspring of body size x surviving to maturity, 
# i.e., theta(x)
x <- 4 * gamma; x  # why are we multiplying by gamma?
sigma5 <- alpha * ((x - gamma)/(1 + beta * (x - gamma))); sigma5

# Number of offspring is
delta / x

# Therefore, finite growth rate is:
(delta / x) * sigma5

## Write a function to calculate sigma for a given body size x
f_sigma <- function(x){
  ifelse(x <= gamma, 0, alpha * ((x - gamma)/(1 + beta * (x - gamma))))
}

# Create a vector of body sizes
size_vec <- seq(gamma / 100, 10 * gamma, by = gamma / 100)
head(size_vec)

# Compute sigma using function
sigma_vec <- f_sigma(size_vec)
tail(sigma_vec)

# Plot survival function
plot(sigma_vec ~ size_vec, type = "l", cex = 0.5, lwd = 3, 
     xlim = c(0, max(size_vec)), 
     xlab = "Body size (g)", 
     ylab = "Annual survival")
abline(v = gamma, lty = 3)

## Compute average # of offspring and plot on semi-log scale
f_offspring_no <- function(x) delta / x
offspring_vec <- f_offspring_no(size_vec)
head(offspring_vec)
plot(offspring_vec ~ size, type = "l", cex = 0.5, log = 'y')

## Finally, compute the finite growth rate lambda by multiplying the two functions
lambda_vec <- offspring_vec * sigma_vec

## Plot
plot(lambda_vec ~ size_vec, type = "l", cex = 0.5)
abline(h = 1, lty = 3)

# Find the value of x which maximizes lambda
size_vec[which(lambda_vec == max(lambda_vec))]
