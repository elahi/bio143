# Spatial Variation

# In this R script we will be exploring growth dynamics when there is spatial variation

# PART I: 1D growth dynamics

#parameter definitions
nParticles <- 100 # initial number of particles
nCells <- 21 #number of cells
lRate <- 0.1 # probability of a particle moving to the left
rRate <- 0.1 # probability of a particle moving to the right
tMax <- 20 # number of time points

# First, create a vector with the growth value that is supported at each location 
spGrowth <- c(rep(0.5, 5), rep(0.99, 6), rep(0.5, 3), rep(1.5, 7))
# Visualize the growth rate values over space
plot(spGrowth, xlab = "Cell", ylab = "Growth Rate")

# define a data frame to record the dynamics of the particles in the cells (columns) over time (rows)
spGrowthDF <- data.frame(matrix(NA,nrow = tMax, ncol = nCells))
spGrowthDF[1,] <- 0 #initialize all cells at zero
spGrowthDF[1,ceiling(nCells/2)] <- nParticles #put all particles in the middle

# helper functions for moving particles
nMoving <- function(i) {
  rbinom(1, i, lRate)
}

nMovingLeft <- function(i) {
  rbinom(1, i, lRate/(lRate + rRate))
}

# Do the same simulation as in the previous script, but this time growth depends on space
for (t in 2:tMax) {
  growth_vals <- ceiling(growthDF[t-1,]*spGrowth)
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

plot_dist <- function(DF, rowN){
  points(as.numeric(DF[rowN,]), pch = 16)
  lines(as.numeric(DF[rowN,]), type = "l")
}

par(mar=c(5, 4, 4, 6) + 0.1)


plot(spGrowth, xlab="", ylab="", 
     axes=FALSE, type="b", col="red")
mtext("Growth Rate",side=4,col="red",line=4) 
axis(4, ylim=c(0,2), col="red",col.axis="red",las=1)

axis(1,pretty(range(cVec),21))
mtext("Location",side=1,col="black",line=2.5) 

par(new=TRUE)

plot(as.numeric(growthDF[1,]), pch = 16, xlab = "", ylab = "", axes=FALSE)
lines(as.numeric(growthDF[1,]), type = "l")
mtext("Number of Individuals",side=2,col="black",line=2.5) 
axis(2, ylim=c(0,150),col="black",las=1)


# plot the change in distribution over time

for (i in 1:20) {
  Sys.sleep(0.5)
  plot_dist(growthDF, i)
}

# Questions to discuss (this is part 3 of homework):
# (i) What is your interpretation of the spatial variation in the model set up? (for instance is the population growing at it's initial condition, will movement increase or decrease growth?)
# (ii) For the initial growth conditions we defined (spGrowth), what do you expect the distribution to look like over time?  
# (iii) What is the probability that any particle or its descendants makes it to the "high growth" region (cell 15)? How does this depend on lRate and rRate? (hint: this is equivalent to starting with nParticles = 1)
# (iv) What happens if the rate of movement (lRate and rRate) are decreased?

# Part 4: Predator-Prey dynamics:

# Part 4: Predator-prey dynamics in space

# helper functions for moving particles
nMoving_prey <- function(i) {
  rbinom(1, i, dPrey)
}

nMovingLeft <- function(i) {
  rbinom(1, i, 0.5)
}

nMoving_pred <- function(i) {
  rbinom(1, i, dPred)
}




# parameter definitions:
nPrey <- 100 # initial number of prey
nPred <- 8
dPrey <- 0.5 # probability of a particle moving to the left
dPred <- 0.3 # probability of a particle moving to the right

# Vector of growth values for prey at each location
spPreyGrowth <- c(rep(0.5, 5), rep(0.99, 6), rep(0.5, 3), rep(1.5, 7))
# Visualize the growth rate values over space
plot(spPreyGrowth, xlab = "Cell", ylab = "Growth Rate")
predRate <- 0.05
predGrowth <- 0.001
predDeath <- 0.02

nCells <- 21 #number of cells
tMax <- 20 # number of time points

#create two different data frames, one for the prey, one for the predator
preyDF <- data.frame(matrix(NA,nrow = tMax, ncol = nCells))
preyDF[1,] <- 0 #initialize all cells at zero
preyDF[1,ceiling(nCells/2)] <- nPrey #put all prey in the middle

predDF <- data.frame(matrix(NA,nrow = tMax, ncol = nCells))
predDF[1,] <- 0 #initialize all cells at zero
predDF[1,ceiling(nCells/2)] <- nPred #put all prey in the middle

for (t in 2:tMax) {
  prey_vals <- ceiling(preyDF[t-1,]*spPreyGrowth) # prey intrinsic growth
  prey_caught <-prey_vals*predDF[t-1,]*predRate #prey death due to predation
  prey_vals <- ceiling(prey_vals - prey_caught) # factor in death due to predation
  
  
  pred_vals <- ceiling(predDF[t-1,]*(1-predDeath)) #intrinsic death rate for predator
  pred_growth <- preyDF[t-1,]*predDF[t-1,]*predGrowth
  pred_vals <- ceiling(pred_vals + pred_growth)
  
  #movement of predator and prey
  nMovingPrey <- sapply(prey_vals, nMoving_prey)
  nMovingLeftPrey <- sapply(nMovingPrey, nMovingLeft)
  nMovingRightPrey <- nMovingPrey - nMovingLeftPrey
  
  nMovingPred <- sapply(pred_vals, nMoving_pred)
  nMovingLeftPred <- sapply(nMovingPred, nMovingLeft)
  nMovingRightPred <- nMovingPred - nMovingLeftPred
  
  #enforce boundaries
  # prey
  nMovingRightPrey[1] <- nMovingRightPrey[1] + nMovingLeftPrey[1]
  nMovingLeftPrey[1] <- 0
  
  nMovingLeftPrey[nCells] <- nMovingRightPrey[nCells] + nMovingLeftPrey[nCells]
  nMovingRightPrey[nCells] <- 0
  
  # predators
  nMovingRightPred[1] <- nMovingRightPred[1] + nMovingLeftPred[1]
  nMovingLeftPred[1] <- 0
  
  nMovingLeftPred[nCells] <- nMovingRightPred[nCells] + nMovingLeftPred[nCells]
  nMovingRightPred[nCells] <- 0
  
  #define the particles at the next time point
  preyDF[t, ] <- prey_vals - nMovingPrey + c(tail(nMovingLeftPrey, -1),0) + c(0, head(nMovingRightPrey, -1))
  predDF[t, ] <- pred_vals - nMovingPred + c(tail(nMovingLeftPred, -1),0) + c(0, head(nMovingRightPred, -1))
}


#plot the distributions over time
cVec <- seq(1:nCells)

plot_dist <- function(DF, rowN){
  points(as.numeric(DF[rowN,]), pch = 16)
  lines(as.numeric(DF[rowN,]), type = "l")
}

par(mar=c(5, 4, 4, 6) + 0.1)


plot(spPreyGrowth, xlab="", ylab="", 
     axes=FALSE, type="b", col="red")
mtext("Growth Rate",side=4,col="red",line=4) 
axis(4, ylim=c(0,2), col="red",col.axis="red",las=1)

axis(1,pretty(range(cVec),21))
mtext("Location",side=1,col="black",line=2.5) 

par(new=TRUE)
plot(as.numeric(preyDF[1,]), pch = 16, col = "blue", xlab = "", ylab = "", axes = FALSE)
lines(as.numeric(preyDF[1,]), type = "l", col = "blue")

points(as.numeric(predDF[1,]), pch = 16, xlab = "", ylab = "")
lines(as.numeric(predDF[1,]), type = "l")

mtext("Number of Individuals",side=2, col="black",line=2.5) 
axis(2, ylim=c(0,150),col="black",las=1)


# plot the change in distribution over time

for (i in 1:20) {
  Sys.sleep(0.5)
  plot_dist(preyDF, i)
  #plot_dist(predDF, i)
}

