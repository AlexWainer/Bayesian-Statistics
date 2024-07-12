# Question 2

set.seed(111)
#### Data Generating Functions ####


generate_lost <- function(grid_size, nsims){

  mu_vec  <- c(grid_size/2, grid_size/2)
  sig_mat <- matrix(c(2, 1, 5, 5), 2,2)
  
  dat <- mvrnorm(nsims, mu_vec, sig_mat)
  dat <- round(abs(dat))
  
  prior <- matrix(rep(0,grid_size^2), grid_size, grid_size)
  for (i in 1:NROW(dat)){
    
    if (dat[i,1] < grid_size & dat[i,2] < grid_size){
      prior[dat[i,1], dat[i,2]] <- prior[dat[i,1], dat[i,2]] + 1
    }
    
  }
  prior <- prior/sum(prior)
  return(prior)
}

# Yi 
generate_fisherman <- function(grid_size){
  
  # Function to generate the true location of the lost fisherman.
  # This should not effect the search decision in any way!! It is unknown
  # to the search crew.
  # Args: 
  #       grid_size: the dimensions of the square search grid
  
  
  mu_vec  <- c(grid_size/2, grid_size/2)
  sig_mat <- matrix(c(2, 1, 5, 5), 2,2)
  
  location  <- round(mvrnorm(1, mu_vec, sig_mat))
  true_grid <- matrix(rep(0, grid_size^2), grid_size, grid_size)
  true_grid[location[1], location[2]] <- 1
  
  return(true_grid)
}

# Function to update theta
update <- function(x, y, theta_prior, detect_pr) {
  # Equation 3
  theta_prior[x, y] <- ((1 - detect_pr[x, y]) * theta_prior[x, y]) / 
    (1 - (detect_pr[x, y] * theta_prior[x, y]))
  # Update other thetas
  for(x_other in 1:nrow(theta_prior)){
    for(y_other in 1:ncol(theta_prior)){
      if(x_other != x && y_other != y){
        # Equation 4
        theta_prior[x_other, y_other] <- theta_prior[x_other, y_other] /
          (1 - (detect_pr[x_other, y_other] * theta_prior[x_other, y_other]))
        
      }
    } 
    
    
  }
  # make sure thetas sum to 1
  theta_prior <- theta_prior / sum(theta_prior)
  
  return(theta_prior)
  
}


#### Simulation ####

search_size <- 20
unifs <- runif(search_size^2, min = 0.6, max = 0.9)
detect_pr <- matrix(unifs, ncol = search_size)



fisherman_location <- generate_fisherman(20)
# get coordinates of fishermans location
fisherman_jk <- which(fisherman_location == max(fisherman_location), arr.ind = T)
fisherman_j <- fisherman_jk[1]
fisherman_k <- fisherman_jk[2]

theta_prior <- generate_lost(20, 1000)
theta_prior_1 <- theta_prior

# Store values
x_searched <- numeric(48)
y_searched <- numeric(48)
fisherman_prob <- numeric(48)

for(i in 1:48){
  i_search <- which(theta_prior == max(theta_prior), arr.ind = T)
  x_coord <- i_search[1]
  y_coord <- i_search[2]
  
  x_searched[i] <- x_coord
  y_searched[i] <- y_coord
  fisherman_prob[i] <- theta_prior[fisherman_j, fisherman_k]
  
  if(fisherman_jk[1] == x_coord && fisherman_jk[2] == y_coord){
    # binomial 1 = Bernoulli
    pr_saved <- rbinom(1, 1, detect_pr[x_coord, y_coord])
    
  } else {
    pr_saved <-0
  }
  
  if(pr_saved == 1){
    
    break
  }
  # Update theta based on posterior evidence
  theta_prior <- update(x_coord, y_coord, theta_prior, detect_pr)
}

theta_prior_final <- theta_prior

data_for_plot <- as.data.frame(as.table(as.matrix(theta_prior_1)))
names(data_for_plot) <- c("X", "Y", "Value")

# Create heatmap
heatmap1 <- ggplot(data_for_plot, aes(x = X, y = Y, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "Initial Probabilities", x = "Column", y = "Row", fill = "Probability") +
  coord_fixed()

# Prepare data 
data_for_plot_2 <- as.data.frame(as.table(as.matrix(theta_prior_final)))
names(data_for_plot_2) <- c("X", "Y", "Value")

heatmap2 <- ggplot(data_for_plot_2, aes(x = X, y = Y, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "Final Probabilities", x = "Column", y = "Row", fill = "Probability") +
  coord_fixed()

grid.arrange(heatmap1, heatmap2, ncol = 2)

plot(fisherman_prob, xlab = "Hour", ylab = "Probability")