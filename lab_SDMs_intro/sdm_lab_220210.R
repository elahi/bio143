################################################################################
##' @title Lab: species distribution modeling with GAMs
##' @author Robin Elahi
##' @contact elahi@stanford.edu
##' @date 2022-02-10
##' @log 
################################################################################

##' Exercise from the Stanford course Quantitative Methods for Marine Ecology and Conservation (BIOHOPK143), 
##' taught by Giulio DeLeo, Maurice Goodman, and Richard Grewelle. 
##' Exercise by Goodman

##### PACKAGES, DATA #####
library(tidyverse)
library(here)
library(mgcv)
library(dismo)

theme_set(theme_bw(base_size = 10) + 
            theme(panel.grid = element_blank(), 
                  strip.background = element_blank())
)

#  Load functions and data
source(here("lab_SDMs_intro", "plot_functions.R"))
d <- read.csv(here("lab_SDMs_intro", "EBS_survey_data.csv"))
head(d)

## Smooth terms: temperature
m_c <- gam(flounder ~ s(bottom_temp, k = 5), 
           data = d, family = binomial(link = "logit"))
summary(m_c)

## Spatio-temporal terms
m_s <- gam(flounder ~ s(E_km, N_km, k = 100), 
           data = d, family = binomial(link = "logit"))

m_s_t <- gam(flounder ~ s(E_km, N_km, k = 100) + s(year, k = 10), 
             data = d, family = binomial(link = "logit"))

m_st <- gam(flounder ~ s(E_km, N_km, k = 100) + s(year, k = 10) +
              ti(E_km, N_km, year, d = c(2,1)), 
            data = d, family = binomial(link = "logit"))

## Spatio-temporal plus temperature 
m_st_c <- gam(flounder ~ s(E_km, N_km, k = 100) + s(year, k = 10) +
                          ti(E_km, N_km, year, d = c(2,1)) + s(bottom_temp, k = 5), 
                        data = d, family = binomial(link = "logit"))

## Compare models
AIC(m_c, m_s, m_s_t, m_st, m_st_c) # use m_st_c for visualization

## Plotting
# scale = 0 allows different y axes for models with multiple smooths
plot(m_c, scale = 0)
plot(m_s, scale = 0)
plot(m_s_t, scale = 0)
plot(m_st, scale = 0)
plot(m_st_c, scale = 0)

## Significance of model terms
anova(m_c)
anova(m_st_c)

##### Characterizing model accuracy #####

fitted <- as.numeric(predict(m_st, type = "response"))
fitted_p <- fitted[d$pollock == 1]
fitted_a <- fitted[d$pollock == 0]

ROC <- dismo::evaluate(fitted_p, fitted_a)
plot(ROC, "ROC")

#### Predicting on new data ####

extrap_grid <- read.csv(here("lab_SDMs_intro", "extrapolation_grid.csv"))

# Using the best model
extrap_grid$fitted <- predict(m_st_c, newdata = extrap_grid, type = "response")

map_EBS_grid <- function(E_km, N_km, fill, facet) {
  
  require(tidyverse)
  
  # Bind data into data frame
  plot_data <- data.frame(E_km = E_km, N_km = N_km, fill = fill, facet = facet)

  ### Plot EBS map
  EBS_map <- plot_data %>% 
    ggplot(aes(E_km, N_km, fill = fill)) + 
    scale_fill_viridis_c(option = "magma") +
    facet_wrap(~facet, ncol = 6) + 
    coord_cartesian(xlim = c(0, 1400), ylim = c(6000, 7000)) +
    geom_raster()
  
  ### Clean up and return map
  theme_map(EBS_map)
  
}

with(extrap_grid, map_EBS_grid(E_km, N_km, fill = fitted, facet = year)) + 
  labs(fill = "Probability of occurrence")

#### Center of gravity ####

N_obs <- d %>% 
  group_by(year) %>% 
  summarise(mean_N = weighted.mean(N_km, w = flounder))

N_fit <- extrap_grid %>% 
  group_by(year) %>% 
  summarise(mean_N = weighted.mean(N_km, w = fitted))

N_obs %>% 
  ggplot(aes(year, mean_N)) + 
  geom_point() + 
  geom_line(data = N_fit)

##### Presentation #####

# A slide listing the 5 models you fit (see the “model selection” section) and their AIC scores. Indicate which model fits best according to AIC, and provide the ROC curve and AUC for this model.

## Model fit
# df       AIC
# m_c      4.688944 14000.298
# m_s     61.085579  7946.941
# m_s_t   74.741584  6098.348
# m_st   126.819111  5860.352
# m_st_c 126.547699  4715.488

fitted <- as.numeric(predict(m_st_c, type = "response"))
fitted_p <- fitted[d$flounder == 1]
fitted_a <- fitted[d$flounder == 0]

ROC <- dismo::evaluate(fitted_p, fitted_a)
plot(ROC, "ROC")

ROC@auc
dismo::threshold(ROC)$spec_sens

# A plot of the temperature smoother from the model with just temperature, and from the model with temperature and spatial/temporal/spatio-temporal terms. What did you expect these curves to look like? Did including spatio-temporal covariates change the shape of the temperature curve?

plot(m_c, scale = 0)
plot(m_s, scale = 0)
plot(m_st_c, scale = 0)

## Logit to prob

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

logit2prob(5)
logit2prob(0)
logit2prob(-5)

d %>% 
  group_by(year) %>%
  summarise(mean = mean(bottom_temp)) %>% 
  ggplot(aes(year, mean)) + 
  geom_point() + 
  geom_smooth() + 
  labs(x = "Year", y = "Mean temperature (C)")
  
