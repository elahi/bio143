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

##### ANALYTICAL SOLUTION #####

# See notes for analytical solution, here we substitute our values:
x_opt = gamma + sqrt(gamma / beta)
x_opt

##### USING R FOR SYMBOLIC CALCULATION #####

# Translate function into an R expression
lambda_expr <- expression(alpha * ((x - gamma)/(1 + beta * (x - gamma))) * delta / x)
lambda_expr

# Ask R to compute derivative
lambda_der <- D(lambda_expr, "x")
lambda_der

# Write function to evaluate the derivative
f_lambda_der <- function(x) eval(lambda_der)

# Is the derivative of the function at x_opt zero?
f_lambda_der(x_opt) # 

# Compute and plot the derivative of lambda as a function of body size x
lambda_der_vec <- sapply(size_vec, f_lambda_der)
plot(lambda_vec ~ size_vec, cex = 0.5, col = "blue")
abline(v = c(gamma, 20 * gamma), h = 0, lty = 3)
abline(v = x_opt, col = "red")

##### ALL THREE FUNCTIONS ON THE SAME PLOT #####

par(mfrow =c(1,1), mar=c(6,4,1,8)) # set margin size to accommodate two y axes on right 
plot(sigma_vec ~ size_vec, type="l", cex=0.5, col="green", lwd=2, xlim=c(0,10*gamma))
leg.txt <- c("Sigma", "Offspring number", "Lambda")
legend(0.0004, 0.000004, leg.txt, col=c("green", "red", "blue"), lty=1, lwd=2, bty="n")

par(new=T) #par adds new data to prior plot without overwriting it 
plot(offspring_vec ~ size_vec, axes = F, cex=0.5, col="red", ylab="", lwd=2, type="l", xlim=c(0,10*gamma))
axis(4) #position 4 is the right side of the quartz window
mtext("Offspring number", side = 4, line=-2) # 'line' adjusts position of axis label
par(new=T) #par adds new data to prior plot without overwriting it 
plot(lambda_vec ~ size_vec, axes = F, cex=0.5, col="blue", type = 'l', ylab="", lwd=2, xlim=c(0,10*gamma))
axis(4, line=4) #position 4 is the right side of the quartz window
mtext("Lambda", side = 4, line=6)
abline(v=x_opt, col="black", lty=2)
