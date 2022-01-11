################################################################################
##' @title Bay anchovy example
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

##### DEFINING PARAMETERS #####

# defining vital rate parameters
#  survival
sigma0 <- 10^-4
sigma1 <- 0.3
sigma2 <- 0.6
sigma3 <- 0.4

# fecundity
f1 <- 5000
f2 <- 10000
f3 <- 50000
f4 <- 100000

# Let's create the Leslie matrix, four rows by four columns
# initially populated by zeros

# we first set the number of classes
noc <- 4

M  <- matrix(0, noc, noc) 

# now we need to populate the first row and sub-diagonal positions with values
M[1,1] <- .5*sigma0*f1 # note that "0.5" is to account for 1:1 sex ratio
M[2,1] <- sigma1
M[1,2] <- .5*sigma0*f2
M[3,2] <- sigma2
M[1,3] <- .5*sigma0*f3
M[4,3] <- sigma3
M[1,4] <- .5*sigma0*f4

# check whether this is the matrix we wanted
M

# defining new variables for our model

# initial condition (arbitrary in this exercise)
num_age_1  <- 0 
num_age_2  <- 50
num_age_3  <- 10
num_age_4  <- 0


num_age_1 <- 125
num_age_2 <- 5
num_age_3 <- 89
num_age_4 <- 5000

n1 <- c(num_age_1,num_age_2,num_age_3,num_age_4); n1 
N1 <- sum(n1); N1 # total population at time zero
w1 <- n1/N1; w1   # initial fraction of individuals in each age class

tmax <- 20        # total number of years to run the model

#  create variables to store the result of simulations
nt <- matrix(0, noc, tmax) # population in each age class
Ntot <- numeric(tmax)      # total population size
wt <- matrix(0, noc, tmax) # proportion of individuals in each age class
lamt <- numeric(tmax)      # year-to-year finite growth rate

# set initial conditions
nt[, 1]  <- c(num_age_1, num_age_2, num_age_3, num_age_4)  
Ntot[1] <- sum(nt[, 1])  
wt[, 1]  <- nt[, 1]/Ntot[1]  

##### COMPUTE POPULATION GROWTH OVER TIME #####

# for loop to compute population growth and structure over time
for (ti in 1:(tmax - 1)) {
  
  nt[, ti + 1] <- M %*% nt[, ti]              # use the recursive equation
  Ntot[ti + 1] <- sum(nt[, ti + 1])           # total population size 
  wt[, ti + 1] <- nt[, ti + 1] / Ntot[ti + 1] # population structure 
  lamt[ti + 1] <- Ntot[ti + 1] / Ntot[ti]     # compute year-to-year finite growth rate
  
}

# create an independent variable, years - for plotting
years <- 1:tmax 

# note here we have to transpose our nt vector so the graph will work

# Plot population size, arithmetic scale
matplot(years,t(rbind(Ntot, nt)), type="l", lwd = 2, 
        xlab="Year", ylab="Number of individuals in age class",
        main="Arithmetic scale") 
legend(x=1, y=0.95*max(Ntot), legend=c("Ntot", "n1", "n2", "n3", "n4"), 
       col=c(1,2,3,4,5), lty = c(1,2,3,4,5), lwd = 2)

# Plot population size, logarithmic scale
matplot(years,t(rbind(Ntot, nt)), lwd = 2, type="l", 
        xlab="Year", ylab="Number of individuals in age class", 
        main="Logarithmic scale", log = 'y')
legend(x=1, y=5*max(Ntot), legend=c("Ntot", "n1", "n2", "n3", "n4"), 
       col=c(1,2,3,4,5), lty = c(1,2,3,4,5), lwd = 2)

# Plot age distribution over time
matplot(years,t(wt),type="l", lwd = 2,
        xlab="Year",ylab="Fraction of individuals in age class",
        main="Fractions in all stages")
legend(x=16, y=0.6, legend=c("n1", "n2", "n3", "n4"), 
       col=c(1,2,3,4), lty = c(1,2,3,4), lwd = 2)

# Plot population growth rate
matplot(years,lamt,type="l", lwd = 2,xlab="Year",ylab="Lambda",
        main="Population growth rate")
abline(h=1, lty = 2) # horizontal line for y=1

#par(sf)

sf <- par(mfrow=c(3, 3)) # let's make a 3x3 panel plot 

barplot(wt[,1], main='age distribution', names.arg = c(1:4))
barplot(wt[,2], main='age distribution', names.arg = c(1:4))
barplot(wt[,3], main='age distribution', names.arg = c(1:4))
barplot(wt[,4], main='age distribution', names.arg = c(1:4))
barplot(wt[,5], main='age distribution', names.arg = c(1:4))
barplot(wt[,6], main='age distribution', names.arg = c(1:4))
barplot(wt[,tmax-2], main='age distribution', names.arg = c(1:4))
barplot(wt[,tmax-1], main='age distribution', names.arg = c(1:4))
barplot(wt[,tmax], main='age distribution', names.arg = c(1:4))

par(sf)

# now, let's use R function to compute eigenvalues and eigenvectors
Ei_Gen <- eigen(M)  
eigenvect <- Ei_Gen$vectors # extract eigenvectors
eigenval <- Ei_Gen$values # extract eigenvalues

# the first element is the *maximum eigenvalue* and correspond to the aymptotic
# population growth rate
Re(eigenval[1]) 

# the eigenvector [,1] corresponds to the maximum eigenvalue 
# and provides the long-term, relative population size structure once we  
#  normalized to 1
stable_stage_distribution  <- Mod(eigenvect[,1])/sum(Mod(eigenvect[,1])) 
stable_stage_distribution
sum(stable_stage_distribution) # check for consistency...

##### SHORTER VERSION PLUS ELASTICITY #####

# As above but shorter, plus elasticity and more
vM <-eigen(M); lambda <-vM$values[1]; Re(lambda) # finite growth rate (max eigenvalue)
ps <-vM$vectors[,1]; ps <-ps/sum(ps); Re(ps) # population structure 
rv <-eigen(t(M))$vectors[,1]; rv <-rv/rv[1]; Re(rv) # reproductive value 
se <-outer(rv,ps)/sum(rv*ps); Re(se) # sensitivity
el <-Re((se*M)/lambda); Re(el)*100 # elasticity [%]

# let's map the transition matrix
image(x = 1:noc,y = 1:noc, t(M), main='transiton matrix',  ylim=c(noc, 0.5))

# as fecundity is way larger than survival, to improve visual
# I take the square root of each element of the transition matrix 
image(x = 1:noc,y = 1:noc, t(M^0.5), main='transiton matrix (square root)',  ylim=c(noc, 0.5))

# now we can plot elasticity
image(x = 1:noc,y = 1:noc, t(el), main='elasticity',  ylim=c(noc, 0.5))

##### ASSIGNMENT #####

# Change initial conditions and describe how long the population takes to reach the stable age distribution 

el_original <- el

el_modified <- el
