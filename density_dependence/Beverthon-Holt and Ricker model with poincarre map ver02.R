####################################################################
# BIOHOPK-81, BIOHOPK143-243 by Giulio De Leo
# 
# Density-dependent dynamics for population with descrete generations
# 
# Objectives
#   - define Beverton-Hold and Ricker density-dependent models
#   - analyze their dynamics in phase diagram (state-space) 
#     and as a function of time
#   - understand the stability or lack thereof
#   - understand the concept of deterministic chaos and 
#     sensitivity to initial conditions
#   - draw your first Poincarré  map
####################################################################

# As usual, first removes all previous material from R's memory
rm(list=ls(all=TRUE))

# Change this to FALSE if you want to plot in RStudio's plot area
plotInNewWindow <- FALSE

####################################################################
# Let's start with
# the BEVERTON HOLT MODEL
####################################################################

K      <-  100   # set here the long-term population carrying capacity
lambda <-  18 # set the per-capita finite growth rate

# compute the density dependent parameter of the BH model accordingly
alpha <- (lambda-1)/K  

#  Define here the function to compute recruitment witthBeverthon Holt model
bh <- function(n){

   lambda*n/(1+alpha*n)

  } 

# Plot the cobweb diagram and the dynamics in time
# Create first a new window external to RStudio if you want to do so
if (plotInNewWindow) windows(width=12, height=6, restoreConsole=TRUE)          

# Plot two panels in the same window,
# one on the left and the other on the right
gf <- par(mfrow=c(1,2))  

# We plot the stock-recruitment function in state-space, i.e. N(t), vs N(t+1),
# with axes intersection at (0,0) 

curve(bh,from=0,to=K*1.5,lwd=4,ylab="N(t+1)",xlab="N(t)", axes = FALSE) 
axis(1, pos=0); axis(2, pos=0)

# draw some stright lines with specific meanig
abline(a=0,b=1) # draw the bisetrix (a stright line with slope equal to 1) 
abline(v=K, lty=3, col="darkgrey") #this is the vertical line at N(t)=K
abline(h=K, lty=3, col="darkgrey") #this is the vertical line at N(t)=K
points(K,K,col=3, bg=2, pch=21, cex=2)

# set the parameters necessary to define: 
# - simulation horizon, etc.
# - initial conditions,
tmax <- 50 # simulation time (number of generations)
n0   <-  1    # initial condition (initial population size)
N    <- numeric(tmax); # initialize vector (by setting aside laptop memory)
N[1] <- n0 # and set here  the initial condition


# this is the recursive equation to compute population size with the BH model 
# at each time step, as a function of population size at the previous time step

for(i in 1:(tmax-1)) {

    N[i+1] <- bh(N[i])

    }  

# Now, finally plot the trajectories in the phase (stage-space) diagram 
lines(N[-length(N)],N[-1],type="s",lwd=2,col="red",ylab="N[t+1]",xlab="N[t]") 

# plot population as a function of time 
plot(N,type="s",lwd=2, col="red", ylab="N(t)", xlab="time") 

# restore original graphic setting
par(gf)

#######################################################################
# now simulate the system for increasing value of the finite growth rate lambda
# go back to (approximately) line 20sh and increase the value of lambda and run the simulations again
#######################################################################

#######################################################################
#######################################################################
# Now, let's look at the dynamics of
# the RICKER  MODEL 
####################################################################

#set parameters for model N(t+1) = lambda*N(t)*exp[-beta*N(t)]

K      <- 100.0         # carrying capacity
lambda <-    1.5         # finite growth rate at low density
betap  <- log(lambda)/K # the density dependent parameter

rf <- function(n,lambda){
   lambda * n * exp(-betap * n) # the Ricker model
} 

n0   <- 10              # this is the initial condition
N    <- numeric(tmax)   # create vector 
N[1] <- n0

for(i in 1:(tmax-1)) {
     N[i+1] <- rf(N[i],lambda)
} 

gf     <- par(mfrow=c(1,2))  
x_max  <-  max(N)

# plot the stock recruiment function
curve(rf(x,lambda), ylim=c(0,x_max*1.15),xlim=c(0,x_max*1.15),
       lwd=3,from=0,to=200,ylab="N(t+1)",xlab="N(t)", axes = FALSE) 
axis(1, pos=0); axis(2, pos=0)
abline(a=0,b=1) # plot the bisetrix

# plot the trajectories in the phase diagram
lines(N[-length(N)],N[-1],type="s",lwd=2,col="red",ylab="N[t+1]",xlab="N[t]")  

# plot population as a function of time
plot(N, type="l", col="red",lwd=2,ylab="N(t)",xlab="t") 
points(N,col="blue",lwd=2)
abline(h=K)

par(gf)

####################################################################
# In order to understand the implications of chaotic systems,
# let's analyze the problem of 
# SENSITIVTY TO INITIAL CONDITIONS
####################################################################

# Let's create a new vector and change  the initial condition 
# by a very small value, when lambda is large (e.g. 18)
N2    <- numeric(tmax)
N2[1] <- n0*1.01

for(i in 1:(tmax-1)) {

  N2[i+1] <- rf(N2[i],lambda)
	
#  if(N2[i+1]==N2[i]){break}
}

gf <- par(mfrow=c(1,1))  

# plot population as a function of time
plot(N,type="l",col="red",lwd=1,ylab="N(t)",xlab="t", lty=3) 
points(N, col="darkgreen",lwd=2, pch=16)
abline(h=K, lty = 3)

points(N2,type="l",col="blue",lwd=1,ylab="N(t)",xlab="t", lty=3)
points(N2,col="brown",lwd=2, pch =16)

# now plot the two time series  on the N-N2 plane to check 
# how close they are
plot(N,N2, type ='p')
points(c(0, max(N, N2)), c(0, max(N, N2)), type ='l', lty = 2)

#######################################################################
# plot the trajectories in the phase diagram
plot(N[-length(N)],N[-1],type="p",lwd=1,col="red",ylab="N[t+1]",xlab="N[t]", pch = 16)  
# plot the trajectories in the phase diagram

lines(N2[-length(N2)],N2[-1],type="p", lwd=1, col="blue",ylab="N[t+1]",xlab="N[t]", pch = 1)  

RV <- rnorm(length(N), mean = K, sd = K/3)


##################################
########## POINCARRE MAP##########
##################################

#define a range of value of the finite growth rate
lambda2 <- seq(1,20,by=0.01)
beta    <- 0.02 # keep beta constant

# throw away the first 100 years of the simulation
burn_in_period <-  100 
buffer  <- numeric(burn_in_period)
tmax    <- 50  # keep the last tmax years


rf <- function(n,lambda){lambda*n*exp(-beta*n)} # the ricker's model

N  <- matrix(0,nrow=length(lambda2),ncol=tmax) # initialize the matrix
n0 <- 10

for(j in 1:length(lambda2)) {

  buffer[1] <-  n0

  for(i in 2:burn_in_period)
    {buffer[i]=rf(buffer[i-1],lambda2[j])}	

# save the last value as initial condition for the next tmax years
  N[j,1] <-  buffer[burn_in_period]  

  for(i in 2:tmax)
    {N[j,i] <- rf(N[j,i-1],lambda2[j])}	
}

if (plotInNewWindow) windows(width=6, height=6, restoreConsole=TRUE)           

par(mfrow=c(1,1))
matplot(lambda2,N[,10:50],pch=".",xlab='lambda',ylab='N', lty=1, col="black")

#######################################################################

############################################
# bonus: POINCARRE MAP of the logistic map #
############################################

#define a range of value of the finite growth rate
lambda2=seq(1,20,by=0.01)
beta=0.02 # keep beta constant

burn_in_period = 100 # throw away the first 100 years of the simulation
buffer = numeric(burn_in_period)
tmax=50  # keep the last tmax years


logmap=function(n,lambda){lambda*n*(1-beta*n)} # the LOGISTIC MAP

N=matrix(0,nrow=length(lambda2),ncol=tmax) # initialize the matrix
n0=10

for(j in 1:length(lambda2)) {

  buffer[1] = n0

  for(i in 2:burn_in_period)
    {buffer[i]=rf(buffer[i-1],lambda2[j])}	

  N[j,1] = buffer[burn_in_period]

  for(i in 2:tmax)
    {N[j,i]=rf(N[j,i-1],lambda2[j])}	
}

if (plotInNewWindow) windows(width=6, height=6, restoreConsole=TRUE)           

matplot(lambda2,N[,10:50],pch=".",xlab='lambda',ylab='N', lty=1, col="black")
