# SEIS Stochastic Simulation using the Gillespie Algorithm
# Running trajectories to compute marginal and joint distributions
rm(list=ls())
# 1. Initialize Parameters
# R=5
N      <- 20
alpha  <- 0.9
beta   <- 0.24    
lambda <- 1.2
n_sims <- 1000
t_max  <- 200

# R=0.533
# N      <- 20
# alpha  <- 1.8
# beta   <- 1.5    
# lambda <- 0.8
# n_sims <- 1000
# t_max  <- 200

# Matrices to store the total time spent in each state across all trajectories
# Rows represent Exposed (0 to N), Columns represent Infectious (0 to N)
joint_time <- matrix(0, nrow = N + 1, ncol = N + 1)
rownames(joint_time) <- 0:N
colnames(joint_time) <- 0:N

set.seed(42)

# 2. Run Trajectories
for (sim in 1:n_sims) {
  # Initial conditions
  e <- 1
  i <- 1
  t <- 0
  
  # Iteration Step
  while ((e + i) > 0 && t < t_max) {
    
    # Calculate the transition rates
    r1 <- i * (N - e - i) * lambda / N
    r2 <- e * alpha
    r3 <- i * beta
    
    R_total <- r1 + r2 + r3
    
    # Failsafe: if total rate is 0, no more events can occur (extinction)
    if (R_total == 0) break
    
    # Generate random numbers
    u1 <- runif(1)
    u2 <- runif(1)
    
    # Determine the time to the next event
    tau <- (1 / R_total) * log(1 / u1)
    
    # Record the time spent in the CURRENT state before the transition
    # (e+1 and i+1 because R indexing starts at 1, but states can be 0)
    joint_time[e + 1, i + 1] <- joint_time[e + 1, i + 1] + tau
    
    # Update time
    t <- t + tau
    
    # Determine which event occurs and update state
    if (u2 < (r1 / R_total)) {
      # Event 1 (New Exposure)
      e <- e + 1
    } else if (u2 < ((r1 + r2) / R_total)) {
      # Event 2 (Progression to Infectious)
      e <- e - 1
      i <- i + 1
    } else {
      # Event 3 (Recovery)
      i <- i - 1
    }
  }
}

# 3. Compute Frequencies (Probabilities)
# Normalize by the total transient time across all simulations to get the distributions
total_time <- sum(joint_time)
joint_dist <- joint_time / total_time

# Marginal distributions (summing over rows and columns)
marginal_E <- rowSums(joint_dist)
marginal_I <- colSums(joint_dist)

# Find the maximum values of E and I actually reached to trim empty plot space
max_idx_E <- max(which(marginal_E > 0))
max_idx_I <- max(which(marginal_I > 0))

# 4. Plotting
# Set up a 1x2 plotting grid so all graphs appear side-by-side
par(mfrow = c(1, 2))

# Plot Marginal Distribution of E
barplot(marginal_E[1:max_idx_E], 
        col = "orange", 
        main = "Marginal Dist. of Exposed (E)",
        xlab = "Number of Exposed", 
        ylab = "Frequency")

# Plot Marginal Distribution of I
barplot(marginal_I[1:max_idx_I], 
        col = "red", 
        main = "Marginal Dist. of Infected (I)",
        xlab = "Number of Infected", 
        ylab = "Frequency")

# Reset plotting parameters
par(mfrow = c(1, 1))
# Plot Joint Distribution of (E, I) as a heatmap
image(x = 0:(max_idx_E - 1), 
      y = 0:(max_idx_I - 1), 
      z = joint_dist[1:max_idx_E, 1:max_idx_I],
      col = heat.colors(50)[50:1], # Invert heat colors so higher density is hotter/darker
      main = "Joint Distribution (E, I)",
      xlab = "Number of Exposed (E)",
      ylab = "Number of Infected (I)")
# Add contour lines on top of the heatmap for easier reading
contour(x = 0:(max_idx_E - 1), 
        y = 0:(max_idx_I - 1), 
        z = joint_dist[1:max_idx_E, 1:max_idx_I], 
        add = TRUE)


