rm(list=ls(all=TRUE)) # removes all previous material from R's memory

df = as.data.frame(read.table("data/Tab3_from_Stolen_and_Barlow_2003_A_MODEL_LIFE_TABLE _FOR_BOTTLENOSE_dolphin3.csv", sep=",", header=T, na.string="-9.99"))
tail(df)


# compute survivorship
df$px <- df$nx/df$nx[1]
age <- 1:nrow(df)

plot(age, df$px, type = 'l', ylim = c(0,1.1))

# mean life expectancy
(e0 <- 1/2 + sum(df$px[-1]))

(e5 <- 1/2 + sum(df$px[-(1:5)]/df$px[5]))


mlex <- function(x)  {
  
  ex <- 1/2 + sum(df$px[- (1:(x+1))]/df$px[(x+1)]) 
  return(ex)
} 

mlex(2)

# apply the function to each age

mle <- sapply(age, mlex)

plot(age, mle, type = 'l')


# compute the rep[roductive number
# 

(Ro <- sum(0.5 * df$px * df$fecx ))


LHS_ELeq <- function(lambda) {sum(0.5 * df$fecx * df$px * lambda^-df$x) }

DTELf <- function(lambda) { LHS_ELeq(lambda)-1 }

# find the root and return it as output of the function
# By default, a function returns the output of its last line.
# Use 'return()' to be more explicit
uniroot(DTELf, c(0.1,2))[1]  # return the root



