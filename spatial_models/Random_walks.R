#Diffusion and random walks
library(tibble)
library(tidyr)
library(ggplot2)

# Part I: Simple Random Walks

# parameter definitions:
nParticles <- 200 # initial number of particles
lRate <- 0.1 # probability of a particle moving to the left; corresponds to D_l 
rRate <- 0.1 # probability of a particle moving to the right; corresponds to D_r
tMax <- 100 # number of time points

# define a data frame to record the dynamics of the particles (columns) over time (rows)
diffusionDF <- data.frame(matrix(NA,nrow = tMax, ncol = nParticles))

# define initial location as 0 
diffusionDF[1,] <- 0

#iterate over for multiple time points
for (t in 2:tMax) {
  mVector <- rbinom(nParticles, 1, lRate + rRate)# vector for whether it's moving or not
  mVector <- mVector * (2 * rbinom(nParticles, 1, lRate/(lRate + rRate)) - 1) # adding directionality 
  diffusionDF[t,] <- diffusionDF[t-1,] + mVector 
}
tVec <- seq(1:tMax)
diffusionDF$Time <- tVec

# plot the trajectories 
diffusionDF_reshaped <- diffusionDF %>% pivot_longer(cols = X1:X100, names_to = "Particle", values_to = "Position")
ggplot(diffusionDF_reshaped, aes(Time, Position, color = as.factor(Particle), group = as.factor(Particle))) + 
  geom_line()

# plot mean distance from the origin as a function of time
mean_dist <- rowMeans(abs(diffusionDF))
ggplot() + aes(tVec, mean_dist) +  geom_point() +
  xlab("Time") + ylab("Mean Distance")

# What function fits the time and distance vectors best? 
# We use nls (non-linear least squares) to fit a power-law relationship
nls(mean_dist ~ b*tVec^a, start = list(a = 1, b = 1))

# 1. Questions to discuss (this is part 1 of homework): 
# (i). Assuming the lRate and rRate are equal (in a symmetric random walk), what is the relationship between the mean distance from the origin and the rate of diffusion (lRate/rRate?)
# (ii). What happens if lRate is not equal to rRate? 
# (ii). What is the scaling relationship between mean distance from the origin and time? How is it estimated in the code above, and can you prove this result analytically? 
# (iii). If we had a two-dimensional lattice and a particle could go in one of eight directions (up, down, left, right, and diagonally), do you think we would have the same scaling relationship as in part (ii)? If not, what would it be? 
# (iv). What would happen if we imposed limits on how far the particle can go? What would the distribution of particles go to over a long time in that case? 

#########################################################
# Class instructions: the first break-out room stops here
#########################################################

# Part II: Growth and Random Walks

# In the previous exercise, the particles simply moved around; they did not "reproduce" or "die". Here, we introduce birth-death processes on the lattice:
# We now modify the previous section to allow for exponential growth. (Note we also change the coding structure here, because the number of particles is not bounded, but the cells are!)

# parameter definitions:
nParticles <- 100 # initial number of particles
nCells <- 21 #number of cells
lRate <- 0.5 # probability of a particle moving to the left
rRate <- 0.5 # probability of a particle moving to the right
tMax <- 20 # number of time points
gRate <- 1.05 #growth rate

# define a data frame to record the dynamics of the particles in the cells (columns) over time (rows)
growthDF <- data.frame(matrix(NA,nrow = tMax, ncol = nCells))
growthDF[1,] <- 0 #initialize all cells at zero
growthDF[1,ceiling(nCells/2)] <- nParticles #put all particles in the middle

nMoving <- function(i) {
  rbinom(1, i, lRate)
}

nMovingLeft <- function(i) {
  rbinom(1, i, lRate/(lRate + rRate))
}

for (t in 2:tMax) {
  growth_vals <- ceiling(growthDF[t-1,]*gRate)
  nMovingT <- sapply(growth_vals, nMoving)
  nMovingLeftT <- sapply(nMovingT, nMovingLeft)
  nMovingRightT <- nMovingT - nMovingLeftT
  
  #enforce boundaries
  nMovingRightT[1] <- nMovingRightT[1] + nMovingLeftT[1]
  nMovingLeftT[1] <- 0
  
  nMovingLeftT[nCells] <- nMovingRightT[nCells] + nMovingLeftT[nCells]
  nMovingRightT[nCells] <- 0
  
  #define the particles at the next time point
  growthDF[t, ] <- growth_vals - nMovingT + c(tail(nMovingLeftT, -1),0) + c(0, head(nMovingRightT, -1))
}

#plot the distributions over time
cVec <- seq(1:nCells)

plot_dist <- function(rowN){
  points(as.numeric(growthDF[rowN,]), pch = 16)
  lines(as.numeric(growthDF[rowN,]), type = "l")
}

plot(as.numeric(growthDF[1,]), pch = 16)
lines(as.numeric(growthDF[1,]), type = "l")
for (i in 1:20) {
  Sys.sleep(0.5)
  plot_dist(i)
}

#plot the mean distance over time
dVec <- abs(cVec - ceiling(nCells/2))
meanDist <- function(row_vals){
  dist_vect <- sum(as.numeric(row_vals)*dVec)/sum(row_vals)
}

mean_dist <- apply(growthDF[,1:nCells], MARGIN = 1, meanDist)
ggplot() + aes(tVec, mean_dist) +  geom_point() +
  xlab("Time") + ylab("Mean Distance")

# What function fits the time and mean distance vectors best? 
# We use nls (non-linear least squares) to fit a power-law relationship
nls(mean_dist ~ b*tVec^a, start = list(a = 1, b=1))



# Questions to discuss (this is part 2 of homework):
# (i). What is the scaling relationship between mean distance from the origin and time now? Can you explain this change intuitively? 
# (ii). What do you expect the density of particles/individuals to look like over a long time? How does this depend on the gRate, lRate and rRate parameters? (You can change some of these parameter values to see what happens)

