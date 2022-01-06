################################################################################
##' @title Sea otters
##' @author Robin Elahi
##' @contact elahi@stanford.edu
##' @date 2022-01-06
##' @log 
################################################################################

##### PACKAGES, DATA #####
library(tidyverse)
library(here)
here()

theme_set(theme_bw(base_size = 10) + 
            theme(panel.grid = element_blank(), 
                  strip.background = element_blank())
)

d <- read.csv(here("sea_otter_growth_rate/data/California_Sea_otters.csv"))
d

##### LIONFISH #####

lambda <- 1.12 # Morris et al. 2011, Bioinvasions
N <- 21.1 # individuals per hectare; Whitfield et al. 2007

## Use for loop to:

# Calculate lionfish density 1 year from now
lambda * N

# Calculate lionfish density 100 years from now

# Vector of population size
n_years <- 101 # 100 years from the start, in this case start = 1
n_vec <- vector(length = n_years)
n_vec[1] <- N

for(i in 2:n_years){
  n_vec[i] <- lambda * n_vec[i - 1]
}

n_vec[101]

lambda^(100) * N

d <- tibble(year = 1:n_years, n = n_vec)

# How many years will it take for the population to double?
d # approximately in year 7 (so 6 years later)
log(2) / log(lambda) # exact answer, 6.1 years

# 10-fold increase 
d %>% print(n = 30) #will take about 21 years

# Plot the results as a function of time (natural scale and semi-log scale)

# Natural scale
d %>% 
  ggplot(aes(x = year, y = n)) + 
  geom_point()

# Semi-log scale
d %>% 
  ggplot(aes(x = year, y = log(n))) + 
  geom_point()

## Question from lab
# Based on what we know about the finite growth rate of the population, how long did the population take to increase from 2 to 22.1ind/ha?

N0 <- 2
tmax <- 50
nt <- numeric(tmax)
nt[1] <- N0

for (t in 1:(tmax - 1)) { 
  nt[t + 1] <- lambda * nt[t]
}

nt
d <- tibble(year = 1:tmax, density = nt)
d %>% 
  ggplot(aes(x = year, y = density)) + 
  geom_step() + 
  geom_point(color = 'red') + 
  geom_line(color = "darkgray", linetype = 3) + 
  labs(x = "Time (years)", y = "Density (ind. / ha)")

# Semi-log scale
d %>% 
  ggplot(aes(x = year, y = density)) + 
  geom_step() + 
  geom_point(color = 'red') + 
  geom_line(color = "darkgray", linetype = 3) + 
  labs(x = "Time (years)", y = "Density (ind. / ha)") + 
  scale_y_log10()

d %>% print(n = 25)

# Create multiplier
mult <- 22.1 / 2
mult

log(mult) / log(lambda)
