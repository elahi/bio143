################################################################################
##' @title Elasmobranch example
##' @author Robin Elahi
##' @contact elahi@stanford.edu
##' @date 2022-01-11
##' @log 
################################################################################

##' Exercise from the Stanford course Quantitative Methods for Marine Ecology and Conservation (BIOHOPK143), 
##' taught by Giulio DeLeo, Maurice Goodman, and Richard Grewelle. 

# Source: Mollet and Cailliet 2002. Mar Freshwater Res. Comparative population
# demography of elasmobranch using life history tables, Leslie matrices and
# stage-based matrix models

##### PACKAGES, DATA #####
library(tidyverse)
library(here)

theme_set(theme_bw(base_size = 10) + 
            theme(panel.grid = element_blank(), 
                  strip.background = element_blank())
)

##### FUNCTIONS #####

# Define a function to compute basic demographic parameters
basic_demographic_parameters <- function(M) {
  
  # finite growth rate (max eigenvalue)
  vM     <- eigen(M); 
  lambda <- Re(vM$values[1]);       # Real part of the first eigenvector
  
  # population structure, i.e., the normalized right eigenvector 
  # associated to the maximum eigenvalue 
  ps  <-  Re(vM$vectors[,1]);   
  ps  <-  ps / sum(ps);             # normalize
  
  # reproductive value, i.e, the left eigenvector associated to the
  # maximum eigenvalue (it is computed as the right eigenvector of the
  # transposed matrix)
  rv  <-  Re(eigen(t(M))$vectors[,1]);
  rv  <-  rv / rv[1];               # standardize relative to 1st class
  
  #Sensitivity and elasticity
  sens <- outer(rv,ps) / sum(rv*ps) # sensitivity
  elas <- (sens*M) / lambda         # elasticity [%]
  
  return(list(lambda=lambda, ps=ps, rv=rv, elast=elas))
  
}

#############################################################################
# let's look first at the full matrix

# create a 10x10 matrix and assign 0 to all the elements
transM <- matrix(nrow=10, ncol=10, byrow=TRUE, data=0)

# set reproductive rate
# Only age 3 and older reproduce
for (j in 3:10) transM[1,j] <- 1.8929

# set survivals
# Fill in the sub-diagonal for ages 1-9
for (j in 1:9) transM[j+1,j] <- 0.631

#check the transition matrix
transM

#plot the transition matrix 
image(1:nrow(transM), 1:ncol(transM), t(transM), ylim =c(ncol(transM)+0.5, 0.5) )


# Call our function
bdp_transM_L <- basic_demographic_parameters(transM); bdp_transM_L 

sum(bdp_transM_L$ps[3:10])

image(1:nrow(bdp_transM_L$elas), 1:ncol(bdp_transM_L$elas), t(Mod(bdp_transM_L$elas)), ylim =c(ncol(transM)+0.5, 0.5))

#############################################################################
# let's now look  at the reduced matrices

reducedM <- matrix(nrow=3, ncol=3, byrow=TRUE, 
                   data<-c(    0,     0, 1.8929,
                               0.631,     0,     0,
                               0,  0.631, 0.6272))

basic_demographic_parameters(reducedM)

############################################################################

reduced2M <- matrix(nrow=2, ncol=2, byrow=TRUE, 
                    data<-c(    0.4104,     1.8929,
                                0.2206,     0.6272))

basic_demographic_parameters(reduced2M)
