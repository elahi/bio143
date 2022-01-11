################################################################################
##' @title Life table example
##' @author Robin Elahi
##' @contact elahi@stanford.edu
##' @date 2022-01-11
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

df <- read_csv(here("age_size", "data", "dolphins_stolen_barlow_2003.csv"))

# Plot trend
plot(age, df$px, type = 'l', ylim = c(0,1.1))

# compute survivorship
df$px <- df$nx / df$nx[1]
age <- 1:nrow(df)

# mean life expectancy
(e0 <- 1/2 + sum(df$px[-1]))
(e5 <- 1/2 + sum(df$px[-(1:5)]/df$px[5]))

mlex <- function(x){
  ex <- 1/2 + sum(df$px[- (1:(x+1))]/df$px[(x+1)]) 
  return(ex)
} 

mlex(2)

# apply the function to each age
mle <- sapply(age, mlex)
plot(age, mle, type = 'l')

# compute the reproductive number
(Ro <- sum(0.5 * df$px * df$fecx ))

LHS_ELeq <- function(lambda) {sum(0.5 * df$fecx * df$px * lambda^df$x) }

DTELf <- function(lambda) { LHS_ELeq(lambda)-1 }

# find the root and return it as output of the function
# By default, a function returns the output of its last line.
# Use 'return()' to be more explicit
uniroot(DTELf, c(0.1,2))[1]  # return the root
