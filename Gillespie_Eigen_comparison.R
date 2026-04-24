# Algorithm 1: Recursive eigenvalue scheme for SEIS Epidemic Model
# Computes the Quasi-Stationary Distribution (v) and dominant eigenvalue (-theta)

# 1. Initialize Parameters
N      <- 20
alpha  <- 0.9
beta   <- 0.24    
lambda <- 1.2
M      <- N + 1 # Number of block levels (e = 0 to N)

# Helper function to get the number of states at each level e
# For e = 0: i goes from 1 to N (Size = N)
# For e > 0: i goes from 0 to N-e (Size = N - e + 1)
get_size <- function(e) {
  if (e == 0) return(N)
  else return(N - e + 1)
}

# 2. Build the Block-Tridiagonal Matrices
Q_diag <- list() # Main diagonal blocks: Q_{e, e}
Q_sub  <- list() # Sub-diagonal blocks: Q_{e, e-1}
Q_sup  <- list() # Super-diagonal blocks: Q_{e-1, e}

for (e in 0:N) {
  k <- e + 1
  size_e <- get_size(e)
  Q_diag[[k]] <- matrix(0, nrow = size_e, ncol = size_e)
  
  # Main Diagonal Blocks: Q_{e, e}
  if (e == 0) {
    for (idx in 1:size_e) {
      i <- idx
      Q_diag[[k]][idx, idx] <- -(i * (N - e - i) * lambda / N + i * beta)
      if (i > 1) Q_diag[[k]][idx, idx - 1] <- i * beta
    }
  } else {
    for (idx in 1:size_e) {
      i <- idx - 1
      Q_diag[[k]][idx, idx] <- -(i * (N - e - i) * lambda / N + e * alpha + i * beta)
      if (i > 0) Q_diag[[k]][idx, idx - 1] <- i * beta
    }
  }
  
  # Super-diagonal Blocks: Q_{e-1, e}
  if (e > 0) {
    size_prev <- get_size(e - 1)
    Q_sup[[k - 1]] <- matrix(0, nrow = size_prev, ncol = size_e)
    
    for (idx_prev in 1:size_prev) {
      i <- if (e - 1 == 0) idx_prev else idx_prev - 1
      idx_curr <- i + 1
      if (idx_curr <= size_e) {
        Q_sup[[k - 1]][idx_prev, idx_curr] <- i * (N - (e - 1) - i) * lambda / N
      }
    }
  }
  
  # Sub-diagonal Blocks: Q_{e, e-1}
  if (e > 0) {
    size_prev <- get_size(e - 1)
    Q_sub[[k - 1]] <- matrix(0, nrow = size_e, ncol = size_prev)
    
    for (idx in 1:size_e) {
      i <- idx - 1
      idx_prev <- if (e == 1) i + 1 else i + 2
      if (idx_prev <= size_prev) {
        Q_sub[[k - 1]][idx, idx_prev] <- e * alpha
      }
    }
  }
}

# =========================================================================
# Algorithm 1: Level-Reduction & Inverse Power Method
# =========================================================================

# Step 1: Choose an initial non-zero guess
v <- list()
total_transient_states <- sum(sapply(0:N, get_size))
for (k in 1:M) {
  v[[k]] <- matrix(1 / total_transient_states, nrow = 1, ncol = get_size(k - 1))
}

# Step 2 & 3: Compute Q_star matrices (Backward Sweep)
Q_star <- list()
Q_star[[M]] <- Q_diag[[M]]

for (i in (M - 1):1) {
  inv_Q_star_next <- solve(-Q_star[[i + 1]])
  Q_star[[i]] <- Q_diag[[i]] + Q_sup[[i]] %*% inv_Q_star_next %*% Q_sub[[i]]
}

# Step 4: Compute R and G matrices
R <- list()
G <- list()
for (i in 1:(M - 1)) {
  inv_Q_star_next <- solve(-Q_star[[i + 1]])
  R[[i]] <- Q_sup[[i]] %*% inv_Q_star_next
  G[[i]] <- inv_Q_star_next %*% Q_sub[[i]]
}

# Step 5: Iteration till convergence
max_iter <- 1000
tol <- 1e-10
u <- list()
v_hat <- list()

for (iter in 1:max_iter) {
  v_prev <- v
  
  # 5.1
  u[[M]] <- v_prev[[M]]
  
  # 5.2
  for (i in (M - 1):1) {
    u[[i]] <- v_prev[[i]] + u[[i + 1]] %*% G[[i]]
  }
  
  # 5.3
  v_hat[[1]] <- u[[1]] %*% solve(Q_star[[1]])
  
  # 5.4
  for (i in 1:(M - 1)) {
    v_hat[[i + 1]] <- u[[i + 1]] %*% solve(Q_star[[i + 1]]) + v_hat[[i]] %*% R[[i]]
  }
  
  # 5.5 Scale \hat{v}^{(n)}
  max_val <- max(sapply(v_hat, function(x) max(abs(x))))
  for (k in 1:M) {
    v[[k]] <- v_hat[[k]] / max_val
  }
  
  # Convergence check
  diff <- max(sapply(1:M, function(k) max(abs(v[[k]] - v_prev[[k]]))))
  if (diff < tol) {
    cat("Inverse Power Method converged in", iter, "iterations.\n")
    break
  }
}

# Step 6: Normalize final QSD and compute eigenvalue
sum_v <- sum(sapply(v, sum))
for (k in 1:M) {
  v[[k]] <- v[[k]] / sum_v
}

theta <- 1 / max_val
cat("Dominant eigenvalue (-theta):", -theta, "\n")

# =========================================================================
# Extract Marginal and Joint Distributions from the computed Exact QSD
# =========================================================================
marginal_E <- rep(0, M)
marginal_I <- rep(0, N + 1)
joint_dist <- matrix(0, nrow = N + 1, ncol = N + 1)

for (e in 0:N) {
  k <- e + 1
  marginal_E[k] <- sum(v[[k]])
  
  size_e <- get_size(e)
  for (idx in 1:size_e) {
    i <- if (e == 0) idx else idx - 1
    
    prob_val <- v[[k]][1, idx]
    marginal_I[i + 1] <- marginal_I[i + 1] + prob_val
    joint_dist[e + 1, i + 1] <- prob_val
  }
}

# =========================================================================
# Run Gillespie Simulation to get Empirical Distributions for Comparison
# =========================================================================
cat("Running Gillespie Simulation for comparison...\n")
n_sims <- 1000
t_max  <- 200
joint_time <- matrix(0, nrow = N + 1, ncol = N + 1)

for (sim in 1:n_sims) {
  e <- 1; i <- 1; t <- 0
  while ((e + i) > 0 && t < t_max) {
    r1 <- i * (N - e - i) * lambda / N
    r2 <- e * alpha
    r3 <- i * beta
    
    R_total <- r1 + r2 + r3
    if (R_total == 0) break
    
    u1 <- runif(1); u2 <- runif(1)
    tau <- (1 / R_total) * log(1 / u1)
    
    joint_time[e + 1, i + 1] <- joint_time[e + 1, i + 1] + tau
    t <- t + tau
    
    if (u2 < (r1 / R_total)) { e <- e + 1 } 
    else if (u2 < ((r1 + r2) / R_total)) { e <- e - 1; i <- i + 1 } 
    else { i <- i - 1 }
  }
}

total_time <- sum(joint_time)
joint_dist_gillespie <- joint_time / total_time
empirical_E <- rowSums(joint_dist_gillespie)
empirical_I <- colSums(joint_dist_gillespie)

# =========================================================================
# Plotting Comparisons
# =========================================================================
# Set up a 1x3 plotting grid to show all three distributions side-by-side
par(mfrow = c(1, 2))

# 1. Comparison of Marginal Exposed (E)
plot(0:N, marginal_E, type = "o", col = "orange", lwd = 2, pch = 16,
     main = "Marginal QSD (Exposed)", xlab = "Number of Exposed", ylab = "Probability")
lines(0:N, empirical_E, type = "o", col = "blue", lwd = 2, pch = 17, lty = 2)
legend("topright", legend = c("Exact", "Gillespie"),
       col = c("orange", "blue"), lwd = 2, pch = c(16, 17), lty = c(1, 2), bty = "n")

# 2. Comparison of Marginal Infected (I)
plot(0:N, marginal_I, type = "o", col = "red", lwd = 2, pch = 16,
     main = "Marginal QSD (Infected)", xlab = "Number of Infected", ylab = "Probability")
lines(0:N, empirical_I, type = "o", col = "blue", lwd = 2, pch = 17, lty = 2)
legend("topleft", legend = c("Exact", "Gillespie"),
       col = c("red", "blue"), lwd = 2, pch = c(16, 17), lty = c(1, 2), bty = "n")

# 3. Exact Joint Distribution as a Heatmap
image(x = 0:N, 
      y = 0:N, 
      z = joint_dist,
      col = heat.colors(50)[50:1], # Invert colors so higher density is red/dark
      main = "Exact Joint QSD (E, I)",
      xlab = "Number of Exposed (E)",
      ylab = "Number of Infected (I)")
contour(x = 0:N, 
        y = 0:N, 
        z = joint_dist, 
        add = TRUE)

# Reset plotting parameters
par(mfrow = c(1, 1))

