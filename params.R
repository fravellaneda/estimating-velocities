b0 <- -1.5 # general intercept
b1 <- 2 # scale parameter for the first Gaussian distribution
b2 <- 2 # scale parameter for the second Gaussian distribution
b3 <- 8 # scale parameter for the third Gaussian distribution

# Mean of the bivariate Gaussian distributions
mu1 <- c(0.2,0.8) 
mu2 <- c(0.4,0.2) 
mu3 <- c(0.8,0.5) 

# Covariances matrices for the bivariate Gaussian distribution
s1 <- matrix(c(0.065,0.030,0.030,0.065), nrow=2) 
s2 <- matrix(c(0.065,-0.030,-0.030,0.065), nrow=2)
s3 <- matrix(c(0.065,0,0,0.065), nrow=2)

# Inverse of covariances matrix
si1 <- solve(s1)
si2 <- solve(s2)
si3 <- solve(s3)

lambda_0 <- 5 # General scale parameter
max.t <- 20 # Number of time stamps
h_t = 1/max.t 
espatial.res <- 42 # Spatial resolution

model <- "220" # Model inlaspacetime