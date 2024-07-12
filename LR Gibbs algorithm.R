library(MASS)
library(cubature)
library(ggplot2)
library(gridExtra)


# Question 1
set.seed(2021)
n = 150 # Number of data points
X.c = data.frame(matrix(rnorm(5*n), ncol=5))
colnames(X.c) = c("X1", "X2", "X3", "X4", "X5")
X = as.matrix(cbind(1, X.c)) # Design matrix
e = matrix(rnorm(n), ncol=1) # Errors
beta.true = matrix(c(1, 0, 10, 0, 2, -3), ncol=1)
Y = X%*%beta.true + e # Observations

gibbs <- function(X, Y, N){
  # priors
  k <- ncol(X)
  n <- nrow(X)
  beta_squiggle <- rep(0,k)
  M <- diag(1, k)
  a <- 1
  b <- 1
  
  # Posterior pre-computes
  XtX <- t(X) %*% X
  beta_hat <- solve(XtX) %*% (t(X) %*% Y)
  
  beta_samples <- matrix(NA, N, k)
  sigma_samples <- numeric(N)
  
  for(i in 1:N){
    
    if(i == 1){
      # use mode 1/(1+1) since mean does not exist
      sigma_2 <- 0.5  
    } else {
      sigma_2 <- sigma_samples[i-1]  # Use the sigma from the previous iteration
    }
    # beta posterior
    sigma_beta <- (solve(M + XtX))
    mu_beta <- sigma_beta %*% ((t(X) %*% Y) + M %*% beta_squiggle)
    
    beta <- mvrnorm(n = 1, mu_beta, sigma_2 * sigma_beta)
    
    # sigma posterior
    
    A2 <- t(Y)%*%(Y) - (t(mu_beta)%*%(M + t(X)%*%(X))%*%(mu_beta))
    
    a_new <- a + n / 2
    b_new <- b + A2 /2
    
    sigma <- 1/rgamma(1, a_new, b_new)
    
    beta_samples[i,] <- beta
    sigma_samples[i] <- sigma
  }
  list(beta = beta_samples, sigma = sigma_samples)
}


## Question c)

results <- gibbs(X,Y,N = 50000)

# delete observations for burn-in
results[[1]] <- results[[1]][5001:50000, ]
results[[2]] <- results[[2]][5001:50000]

par(mfrow = c(2, 3))

for(i in 1:ncol(results[[1]])){
  
  coeff_series_i <- ts(results[[1]][, i])
  
  plot_title <- paste("Trace Plot for Coefficient", i-1)
  
  # Plot the time series with specified main title, x label, y label, and color
  plot.ts(coeff_series_i, col = "navy", xlab = "Iteration", ylab = paste("Beta", i-1) , main = plot_title)
}


par(mfrow = c(1, 1))

par(mfrow = c(3, 3))

sample_averages <- apply(results[[1]], 2, mean)

for(i in 1:ncol(results[[1]])){
  
  beta_i <- results[[1]][,i]
  density_i <- density(beta_i)
  plot(density_i, main = paste("Density of Beta", i-1), 
       xlab = paste("Beta", i-1), ylab = "Density")
  
  
  abline(v = sample_averages[i], col="red", lw = 2)
  
  abline(v = beta.true[i], col = "blue", lw = 2)
  
  cred_interval_i <- quantile(beta_i, probs = c(0.025, 0.975))
  abline(v=cred_interval_i, col="darkgreen", lwd=2)
  
  
}

plot(1, type="n", axes=FALSE, xlab="", ylab="", main="") # Create an empty plot
legend("center", legend=c("Sample Average", "True Value", "95% Credible Interval"),
       col=c("red", "blue", "darkgreen"), lwd=2)

par(mfrow = c(1, 1))