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
AIC(m_c, m_s, m_s_t, m_st, m_st_c)

##

