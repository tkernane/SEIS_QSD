# Algorithm 3: The stationary probability of the finite QBD
rm(list=ls())
N      <- 10      # total population
alpha  <- 1.8     # rate E -> I
beta   <- 1.5    # rate I -> S
lambda <- 0.8     # infection rate

alpha_e  <- function(e) alpha * e
beta_i   <- function(i) beta * i
lambda_ei <- function(e, i) i * (N - e - i) * lambda / N

#----------------------------------------------
# The block Q_{e,e-1}
Q_em1 <- function(e) {
  if (e == 1) {
    return( diag(rep(alpha_e(1), N)) )   # Case e = 1 : Q_{1,0}
  }
  m <- N - e + 1                         # Case e > 1 : Q_{e,e-1}
  Q <- matrix(0, m, m + 1)
  
  for (i in 1:m) {
    Q[i, i + 1] <- alpha_e(e)
  }
  return(Q)
}

#-------------------------
# The block Q_{e,e}
Q_ee <- function(e, beta_fun = beta_i) {
  # =========================
  # Case e = 0 : Q_{0,0}
  # =========================
  if (e == 0) {
    Q <- matrix(0, N, N)
    for (i in 1:N) {
      q0i <- beta_fun(i) + lambda_ei(0, i)
      Q[i, i] <- -q0i
      if (i >= 2) {
        Q[i, i - 1] <- beta_fun(i)
      }
    }
    return(Q)
  }
  
  # =========================
  # Case e >= 1 : Q_{e,e}
  # =========================
  m <- N - e + 1
  Q <- matrix(0, m, m)
  
  # FIXED: Added parentheses around (m - 1) to prevent i from starting at -1
  for (i in 0:(m - 1)) {
    q_ei <- alpha_e(e) + beta_fun(i) + lambda_ei(e, i)
    Q[i+1, i+1] <- -q_ei
    if (i >= 1) {
      Q[i+1, i ] <- beta_fun(i)
    }
  }
  return(Q)
}

#-----------------------
# The block Q_{e,e+1}
Q_ep1<- function(e) {
  if (e == 0) {                            # Case e = 0 : Q_{0,1}
    Q <- matrix(0, N, N)
    for (i in 1:(N - 1)) {
      Q[i, i + 1] <- lambda_ei(0, i)
    }
    return(Q)
  }
  
  m1 <- N - e + 1                         # Case e > 0 : Q_{e,e+1}
  m2 <- N - e
  Q <- matrix(0, m1, m2)
  for (i in 0:(m2 - 1)) {
    Q[i + 1, i + 1] <- lambda_ei(e, i)
  }
  return(Q)
}

#----------------------------------------------------------------------------
# ALGORITHM of GAVER
#----------------------------------------------------------------------------
# Step 1: compute matrices B_e
#.............................................................................
compute_B <- function(N) {
  B <- vector("list", N + 1)
  B[[1]] <- Q_ee(0)
  
  for (e in 1:N) {
    B[[e + 1]] <- Q_ee(e) +
      Q_em1(e) %*% solve(-B[[e]]) %*% Q_ep1(e - 1)
  }
  B
}

#-----------------------------------------------------------------------------
# Step 2: compute P(N)B_N=0 and P(N)e_d(N)=1
#------------------------------------------------------------------------------
stationary_dist <- function(Q) {
  d <- nrow(Q)
  A <- t(Q)
  A[d, ] <- rep(1, d)
  b <- c(rep(0, d - 1), 1)
  
  solve(A, b)
}

#-------------------------------------------------------------------------------
# Steps 3 and 4: decreasing reconstruction + normalization
#-------------------------------------------------------------------------------
Gaver_QBD <- function(N) {
  B <- compute_B(N)
  P <- vector("list", N + 1)
  
  P[[N + 1]] <- stationary_dist(B[[N + 1]])
  delta <- sum(P[[N + 1]])
  
  for (e in (N - 1):0) {
    P[[e + 1]] <- as.numeric(
      P[[e + 2]] %*% Q_em1(e + 1) %*% solve(-B[[e + 1]])
    )
    delta <- delta + sum(P[[e + 1]])
  }
  
  for (e in 0:N) P[[e + 1]] <- P[[e + 1]] / delta
  P
}

#----------------------------------------------------------------------------
#=====================================
# Stationary joint distribution
# rows : Exposed E
# Columns : Infected I
#=====================================
joint_dist <- function(P, N) {
  Pi <- matrix(0, nrow = N + 1, ncol = N + 1)
  
  # ---- Niveau e = 0 (I = 1,...,N)
  Pi[1, 2:(N + 1)] <- as.numeric(P[[1]])
  
  # ---- Niveaux e >= 1 (I = 0,...,N-e)
  for (e in 1:N) {
    Pi[e + 1, 1:(N - e + 1)] <- as.numeric(P[[e + 1]])
  }
  
  rownames(Pi) <- paste0("E=", 0:N)
  colnames(Pi) <- paste0("I=", 0:N)
  
  return(Pi)
}

#----------------------------------------------------------------------------


P  <- Gaver_QBD(N)
QBD_joint_dist <- joint_dist(P, N)

marg_E <- rowSums(QBD_joint_dist)
marg_I <- colSums(QBD_joint_dist)

# =========================================================================
# Run Gillespie Simulation to get Empirical Distributions for Comparison
# =========================================================================
cat("Running Gillespie Simulation for comparison...\n")
set.seed(42) # Added a seed for reproducible trajectories 
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
# Set up a 1x2 plotting grid to show all three distributions side-by-side
par(mfrow = c(1, 2))

# 1. Comparison of Marginal Exposed (E)
plot(0:N, marg_E, type = "o", col = "orange", lwd = 2, pch = 16,
     main = "Marginal QSD (Exposed)", xlab = "Number of Exposed", ylab = "Probability")
lines(0:N, empirical_E, type = "o", col = "blue", lwd = 2, pch = 17, lty = 2)
legend("topright", legend = c("QBD", "Gillespie"),
       col = c("orange", "blue"), lwd = 2, pch = c(16, 17), lty = c(1, 2), bty = "n")

# 2. Comparison of Marginal Infected (I)
plot(0:N, marg_I, type = "o", col = "red", lwd = 2, pch = 16,
     main = "Marginal QSD (Infected)", xlab = "Number of Infected", ylab = "Probability")
lines(0:N, empirical_I, type = "o", col = "blue", lwd = 2, pch = 17, lty = 2)
legend("topright", legend = c("QBD", "Gillespie"),
       col = c("red", "blue"), lwd = 2, pch = c(16, 17), lty = c(1, 2), bty = "n")

par(mfrow = c(1, 1))

# 3. QBD Joint Distribution as a Heatmap
image(x = 0:N, 
      y = 0:N, 
      z = QBD_joint_dist,
      col = heat.colors(50)[50:1], # Invert colors so higher density is red/dark
      main = "Joint QSD of the Finite QBD Process (E, I)",
      xlab = "Number of Exposed (E)",
      ylab = "Number of Infected (I)")
contour(x = 0:N, 
        y = 0:N, 
        z = QBD_joint_dist, 
        add = TRUE)
