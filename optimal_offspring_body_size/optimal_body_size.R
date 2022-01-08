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

