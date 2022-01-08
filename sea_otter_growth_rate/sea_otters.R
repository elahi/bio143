################################################################################
##' @title Estimating sea otter population growth rate
##' @author Robin Elahi
##' @contact elahi@stanford.edu
##' @date 2022-01-06
##' @log 
################################################################################

##' Example from Stanford course "Quantitative Methods in Marine Ecology and Conservation"
##' 2022: taught by Giulio DeLeo, Maurice Goodman, Richard Grewelle

##### PACKAGES, DATA #####
library(tidyverse)
library(here)

theme_set(theme_bw(base_size = 10) + 
            theme(panel.grid = element_blank(), 
                  strip.background = element_blank())
)

d <- read_csv(here("sea_otter_growth_rate", "data", "California_Sea_otters.csv"))
d

##### ESTIMATE AVERAGE SEA OTTER POPULATION GROWTH RATE, ASSUMING MALTHUSIAN GROWTH #####

# Plot
d %>% 
  ggplot(aes(x = year, y = abundance)) + 
  geom_step() + 
  geom_point(color = 'red') + 
  geom_line(color = "darkgray", linetype = 3) + 
  labs(x = "Time (years)", y = "Abundance") 

# We would like to estimate lambda
# and add a trendline and confidence interval to the above plot

# Estimate lambda using linear model with log-transformed abundance
d <- d %>% 
  mutate(log_abundance = log(abundance))

m1 <- lm(log_abundance ~ year, data = d)
summary(m1)

# Parameter coefficients
coef(m1)

# Assign slope to an object (igr; instantaneous growth rate)
igr <- coef(m1)[2]
igr

# Take the exponent to get lambda (i.e., finite growth rate)
lambda <- exp(igr)
lambda

# Extract confidence intervals
confint(m1, 'year', level = 0.95)

# Create a new dataframe for predictions and plotting
# Would like to make a prediction every year (this is optional; can also use the original data, d)
d_new <- tibble(year = min(d$year):max(d$year)) 
d_new

# Get predictions from model
d_pred <- predict(m1, newdata = d_new, interval = "confidence") 
d_pred
glimpse(d_pred) # not a dataframe yet

d_pred2 <- d_new %>% 
  mutate(fit = exp(d_pred[, 'fit']), 
         lower = exp(d_pred[, 'lwr']), 
         upper = exp(d_pred[, 'upr']))

d_pred2

# Plot
d_pred2 %>% 
  ggplot(aes(x = year, y = fit)) + 
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) + 
  geom_point(data = d, aes(x = year, y = abundance), color = 'red') + 
  labs(x = "Time (years)", y = "Abundance") 

##### EXERCISE: GRAY SEAL PUPS #####

# Change the input data, rename to match the otter data
d <- read_csv(here("sea_otter_growth_rate", "data", "Grey_seal_puppy_production_Sable_Island_time_trend_data.csv"))
d

d <- d %>% 
  mutate(year = Year)

# Plot
d %>% 
  ggplot(aes(x = year, y = abundance)) + 
  geom_step() + 
  geom_point(color = 'red') + 
  geom_line(color = "darkgray", linetype = 3) + 
  labs(x = "Time (years)", y = "Abundance") 

# We would like to estimate lambda
# and add a trendline and confidence interval to the above plot

# Estimate lambda using linear model with log-transformed abundance
d <- d %>% 
  mutate(log_abundance = log(abundance))

m1 <- lm(log_abundance ~ year, data = d)
summary(m1)

# Parameter coefficients
coef(m1)

# Assign slope to an object (igr; instantaneous growth rate)
igr <- coef(m1)[2]
igr

# Take the exponent to get lambda (i.e., finite growth rate)
lambda <- exp(igr)
lambda

# Extract confidence intervals
confint(m1, 'year', level = 0.95)

# Create a new dataframe for predictions and plotting
# Would like to make a prediction every year (this is optional; can also use the original data, d)
d_new <- tibble(year = min(d$year):max(d$year)) 
d_new

# Get predictions from model
d_pred <- predict(m1, newdata = d_new, interval = "confidence") 
d_pred
glimpse(d_pred) # not a dataframe yet

d_pred2 <- d_new %>% 
  mutate(fit = exp(d_pred[, 'fit']), 
         lower = exp(d_pred[, 'lwr']), 
         upper = exp(d_pred[, 'upr']))

d_pred2

# Plot
d_pred2 %>% 
  ggplot(aes(x = year, y = fit)) + 
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) + 
  geom_point(data = d, aes(x = year, y = abundance), color = 'red') + 
  labs(x = "Time (years)", y = "Abundance") 

##### EXERCISE: LOGGERHEAD TURTLES #####

# Change the input data, rename to match the otter data
d <- read_csv(here("sea_otter_growth_rate", "data", "australian_gbr_Loggerhead_sea_turtle_data.csv"))
d

d <- d %>% 
  mutate(abundance = count)

# Plot
d %>% 
  ggplot(aes(x = year, y = abundance)) + 
  geom_step() + 
  geom_point(color = 'red') + 
  geom_line(color = "darkgray", linetype = 3) + 
  labs(x = "Time (years)", y = "Abundance") 

# We would like to estimate lambda
# and add a trendline and confidence interval to the above plot

# Estimate lambda using linear model with log-transformed abundance
d <- d %>% 
  mutate(log_abundance = log(abundance))

m1 <- lm(log_abundance ~ year, data = d)
summary(m1)

# Parameter coefficients
coef(m1)

# Assign slope to an object (igr; instantaneous growth rate)
igr <- coef(m1)[2]
igr

# Take the exponent to get lambda (i.e., finite growth rate)
lambda <- exp(igr)
lambda

# Extract confidence intervals
confint(m1, 'year', level = 0.95)

# Create a new dataframe for predictions and plotting
# Would like to make a prediction every year (this is optional; can also use the original data, d)
d_new <- tibble(year = min(d$year):max(d$year)) 
d_new

# Get predictions from model
d_pred <- predict(m1, newdata = d_new, interval = "confidence") 
d_pred
glimpse(d_pred) # not a dataframe yet

d_pred2 <- d_new %>% 
  mutate(fit = exp(d_pred[, 'fit']), 
         lower = exp(d_pred[, 'lwr']), 
         upper = exp(d_pred[, 'upr']))

d_pred2

# Plot
d_pred2 %>% 
  ggplot(aes(x = year, y = fit)) + 
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) + 
  geom_point(data = d, aes(x = year, y = abundance), color = 'red') + 
  labs(x = "Time (years)", y = "Abundance") 

