# Density dependent dynamics, the logistic model
# BIOHOPK-81 - BIOHOPK143
#
# removes all previous variables from R's memory
rm(list=ls(all=TRUE)) 

# To simulate the dynamics of logisticgrowth,  we need to define:
# - the model and its model parameters
# - initial condition
# - time horizon
# - time steps at which we want the results of model integration 

#################################################################
# vector of labeled model parameters
mod.pars <- c(r = 0.1, K = 100); mod.pars
r   # try to print out 'r'... 
mod.pars['r']  # extract the value of 'r'

# there is an alternative way that uses the function 'with'
# it seems clumsy, but it is convenient when you have 
# lots of parameters 
with(as.list(mod.pars),print(r)) 


# Set up the function to describe the logistic model
# It is written with a very specific format required by the function
# 'lsoda' in package 'deSolve' to integrate this dynamics model in time
# 
#            lsoda(Init_conditions, Time, model, params)
# the input parameters of the model function are 
# time 't', state variable(s) 'N' and a vector for named parameters
# of the model, i.e., 'params'  

logistic_model  <- function(t, N, params) { 
  
  with(as.list(params), {  # "with" extract named parameters from "params" 
    
    dN.dt <- r * (1 - N / K) * N
    
    return(list(dN.dt))    # 'return' returns the results in a list
    
  })
}

#################################################################
# define the vector of time 
Tmax     <- 100 # time horizon  
TimeStep <- 1 # integration time step
Time     <- seq(0, Tmax, by = TimeStep)  # the corresponding vector


#################################################################
# We will use the function lsoda from package deSolve to integrate the ODE 

if (!require("deSolve")) install.packages("deSolve")
library(deSolve)

Lgstc.output <- lsoda(c(PopSize = 1), Time, logistic_model, mod.pars)

head(Lgstc.output)
tail(Lgstc.output)

plot(PopSize~time, data=Lgstc.output, type='l', lwd=3)

#################################################################
# let's now simulate the system with different initial conditions
# 

Lgm <- lsoda(c(1, 10, 25, 50, 75, 110, 125) , Time, logistic_model, mod.pars)

matplot(x = Lgm[, 1], y = Lgm[, 2:ncol(Lgm)], type = 'l', lwd = 2)

abline(h=mod.pars['K'], col='red', lty=2)
abline(h=mod.pars['K']/2, col='gray',lty=3)

#################################################################
# let's plot the stock recruitment function
attach(as.list(mod.pars))    # provide access to all the parameters

Nv   <- seq(0, K*1.1, by=1)  # create a sequence of density values


dNdt <-  logistic_model(0, Nv, mod.pars)[[1]] # compute production (RHS)

plot(dNdt~Nv, type='l', lwd=2)                # plot it
abline(h=0)
points(0~0, type='p', pch=0, lwd=2, cex=1.5, col = 'darkred')
points(0~K, type='p', pch=16, lwd=2, cex=2, col = 'darkred')
abline(v=K/2, lty=3)

detach(as.list(mod.pars))    # detach the list

#################################################################
generalized_LGM  <- function(t, N, params) { 
  
  with(as.list(params), {  # "with" extract named parameters from "params" 
    # dN.dt <- r *  (1 - N / K)      *  N          # This is the logistic model
    # dN.dt <- r *  ((1 - N/K)^alpha * (N^beta)    # This is the generalized logistic model
    dN.dt <- r * Re(as.complex(1 - N/K)^alpha)* (N^beta)
    return(list(dN.dt))    # 'return' returns the results in a list
  })
}

mod.pars <- c(r = 0.1, K = 100, alpha = 1.2, beta = 1); mod.pars

gLG <- lsoda(c(1, 10, 25, 50, 75, 110, 125) , Time, generalized_LGM, mod.pars)

matplot(x = gLG[, 1], y = gLG[, 2:ncol(gLG)], type = 'l')

abline(h=mod.pars['K'],   col='red', lty=2)
abline(h=mod.pars['K']/2, col='gray',lty=3)


#################################################################
# let's plot the stock recruitment function

mod.pars['alpha']  <- 1; mod.pars['beta']   <- 1

mod.pars

attach(as.list(mod.pars))    # provide access to all the parameters

Nv   <- seq(0, K*1.1, by=1)  # create a sequence of density values

dNdt_LGn1 <-  generalized_LGM(0, Nv, mod.pars)[[1]] # compute production (RHS)

mod.pars['alpha']  <- 1.1
dNdt_LGa1 <-  generalized_LGM(0, Nv, mod.pars)[[1]] # compute production (RHS)

mod.pars['alpha']  <- 1.2
dNdt_LGa2 <-  generalized_LGM(0, Nv, mod.pars)[[1]] # compute production (RHS)

mod.pars['alpha']  <- 0.9
dNdt_LGa3 <-  generalized_LGM(0, Nv, mod.pars)[[1]] # compute production (RHS)

mod.pars['alpha']  <- 0.8
dNdt_LGa4 <-  generalized_LGM(0, Nv, mod.pars)[[1]] # compute production (RHS)

mod.pars['alpha']   <- 1
mod.pars['beta']  <- 1.1
dNdt_LGb1 <-  generalized_LGM(0, Nv, mod.pars)[[1]] # compute production (RHS)

mod.pars['beta']  <- 1.2
dNdt_LGb2 <-  generalized_LGM(0, Nv, mod.pars)[[1]] # compute production (RHS)

mod.pars['beta']  <- 0.9
dNdt_LGb3 <-  generalized_LGM(0, Nv, mod.pars)[[1]] # compute production (RHS)

mod.pars['beta']  <- 0.8
dNdt_LGb4 <-  generalized_LGM(0, Nv, mod.pars)[[1]] # compute production (RHS)

y.max <-max(dNdt_LGn1, dNdt_LGa1, dNdt_LGa2, dNdt_LGa3, dNdt_LGa4, 
                       dNdt_LGb1, dNdt_LGb2, dNdt_LGb3, dNdt_LGb4)



plot(dNdt_LGn1~Nv, type='l', lwd=2, ylim = c(0,y.max))                # plot it
abline(h=0)
points(0~0, type='p', pch=0, lwd=2, cex=1.5, col = 'darkred')
points(0~K, type='p', pch=16, lwd=2, cex=2, col = 'darkred')
abline(v=K/2, lty=3)

points(dNdt_LGa1~Nv, type = 'l', lty = 2, col = 'red')
points(dNdt_LGa2~Nv, type = 'l', lty = 2, col = 'red')
points(dNdt_LGa3~Nv, type = 'l', lty = 1, col = 'red')
points(dNdt_LGa4~Nv, type = 'l', lty = 1, col = 'red')

points(dNdt_LGb1~Nv, type = 'l', lty = 2, col = 'green')
points(dNdt_LGb2~Nv, type = 'l', lty = 2, col = 'green')
points(dNdt_LGb3~Nv, type = 'l', lty = 1, col = 'green')
points(dNdt_LGb4~Nv, type = 'l', lty = 1, col = 'green')


detach(as.list(mod.pars))    # detach the list

#################################################################


#################################################################
# Allee effects
#################################################################
# vector of model parameters
Amod.pars    <- c(r = 0.05, K = 100, A=25); mod.pars

Allee_model  <- function(t, N, params) { #back to classic logistic model  
    with(as.list(params), {  # "with" extract named parameters from "params" 
    dN.dt <- r * (N / A - 1) * (1 - N / K) * N
    return(list(dN.dt))    # 'return' returns the results in a list
  })
}


#################################################################
# draw the stock-recruitment function
# 
attach(as.list(Amod.pars))
Nv   <- seq(0,K*1.1,by=1)
dNdt <- Allee_model(0,Nv,mod.pars)[[1]]
plot(dNdt~Nv, type='l', lwd=2)
abline(h=0)
points(0~0, type='p', pch=16, lwd=2, cex=2,   col = 'darkred')
points(0~A, type='p', pch=0,  lwd=2, cex=1.5, col = 'darkred')
points(0~K, type='p', pch=16, lwd=2, cex=2,   col = 'darkred')
abline(v=A, lty=3)
detach(as.list(Amod.pars))
#################################################################

#################################################################
# let's now simulate the system with different initial conditions
# 

Lgm <- lsoda(c(1, 10, 20, 30, 50, 70, 120) , Time, Allee_model, Amod.pars)

matplot(x = Lgm[, 1], y = Lgm[, 2:ncol(Lgm)], type = 'l', lwd= 2)

abline(h = Amod.pars['K'], col='black', lty=2)
abline(h = Amod.pars['K']/2, col='gray',lty=3)
abline(h = Amod.pars['A'], col='red', lty=5)

