library(deSolve)

library(parallel)

library(doParallel)

library(foreach)

# Model name and date/time

options(error = traceback)

MODEL <- "simpleAOP"

DATETIME <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Create directories

pathToOutput <- paste0("/home/ditilliof/Phd material/qAOP_framework/feedback qAOP/",
                       MODEL, "/")

if (!dir.exists(pathToOutput)) {
  
  dir.create(pathToOutput)
  
}

# Define the model

AOPP <- function(t, ic, k, d1, s2, d2, eps) {
  
  KE1 <- ic[1]
  
  KE2 <- ic[2]
  
  dKE1 <- -d1 * KE1 +  k * KE2 + eps
  
  dKE2 <- s2 * (KE1^n) / (1 + (KE1^n)) - d2 * KE2 + eps
  
  return(c(dKE1, dKE2))
  
}

# Print message

cat("The model is defined\n")

# Set parameters

n <- 3

t <- c(0, 10000)

# Define function to solve the model

solveAOP <- function(k, d1, s2, d2, eps) {
  
  result <- ode(y = c(0.2, 0.2), times = t, func = AOPP, parms = c(k, d1, s2, d2, eps))
  
  ic <- tail(result, 1)
  
  result <- deSolve::ode(y = ic, times = t, func = AOPP, parms = c(k, d1, s2, d2, eps))
  
  return(result[2, 2])
  
}

# Print message

cat("Solve function is defined\n")

# Parameter sets

ress <- 5

k_set <- seq(0.10, 0.3, length.out = ress)

d1_set <- seq(0.1, 2, length.out = ress)

s2_set <- seq(3, 5, length.out = ress)

d2_set <- seq(0.1, 2, length.out = ress)

eps_set <- seq(0.1, 0.3, length.out = ress)

k <- 1

# Print message

cat("Sets are defined\n")

# Generate items

items <- expand.grid(k = k_set, d1 = d1_set, s2 = s2_set, d2 = d2_set, eps = eps_set)

# Print message

cat("Items are defined\n")

# Parallel computation

# Parallel computation
# cl <- makeCluster(detectCores())
# clusterExport(cl, c("AOPP", "ode", "t", "solveAOP", "deSolve::ode"))  # Exporting necessary functions
# KE2_vec <- parSapply(cl, items, function(params) solveAOP(params[1], params[2], params[3], params[4], params[5]))
# stopCluster(cl)

# Initialize cluster
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# Define function to solve the model
solveAOP <- function(params) {
  # Your solveAOP implementation here
}

# Perform parallel computation
KE2_vec <- foreach(params = items, .combine = c) %dopar% {
  solveAOP(params)
}

# Stop cluster
stopCluster(cl)


# Normalize results

KE2_max <- max(KE2_vec)

AO_sc <- KE2_vec / KE2_max

# Print message

cat("The time of simulation for the parallel computation is", Sys.time() - Sys.time(), "seconds\n")

# Plotting

pdf(paste0(pathToOutput, DATETIME, "_FDT_AOdistribution.pdf"))

hist(AO_sc, breaks = 50, main = "AO score distribution", xlab = "AO score", col = "lightblue", probability = TRUE)

dev.off()
