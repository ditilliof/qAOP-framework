library(deSolve)
library(tidyverse)

# Resolve conflicts
conflict_prefer("filter", "dplyr")

# Create a timestamped folder
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_folder <- paste0("~/Phd material/qAOP_framework/model_comparison/", timestamp)
dir.create(output_folder, recursive = TRUE)

# Define initial states and parameters
inistateB = c(MIE = 0.5, KE1 = 0, KE2 = 0)
parsB = c(tau = 0.05, eps = 0.2, d1 = 0.5, k_12 = 0.3, d2 = 0.4)

# Generate data from the model by adding noise to the solutions of the system
finish <- seq(1, 100, by = 10)
outB <- as.data.frame(ode(inistateB, finish, func = qAOPB, parms = parsB))

nrep <- 5  # number of fake replicates
std <- 0.1  # standard deviation to add noise to the solution

fakedataB <- data.frame(time = numeric(), MIE = numeric(), KE1 = numeric(), KE2 = numeric(), REPLICATE = numeric())

for (i in 1:length(finish)) {
  for (j in 1:nrep) {
    fakedataB <- rbind(fakedataB, data.frame(time = finish[i], MIE = outB$MIE[i] + rnorm(1, mean = 0, sd = std),
                                             KE1 = outB$KE1[i] + rnorm(1, mean = 0, sd = std), 
                                             KE2 = outB$KE2[i] + rnorm(1, mean = 0, sd = std), REPLICATE = j))
  }
}
write.csv(fakedataB, paste0(output_folder, "/fakedata_qAOPB.csv"))

# Initial state and parameters for model A
inistateA = c(MIE = 0.5, KE1 = 0, KE2 = 0)
parsA = c(tau = 0.05, eps = 0.2, d1 = 0.5, k_12 = 0.3, d2 = 0.4, k_2 = 0.05)
outA = as.data.frame(ode(inistateA, finish, func = qAOPA, parms = parsA))

fakedataA <- data.frame(time = numeric(), MIE = numeric(), KE1 = numeric(), KE2 = numeric(), REPLICATE = numeric())

for (i in 1:length(finish)) {
  for (j in 1:nrep) {
    fakedataA <- rbind(fakedataA, data.frame(time = finish[i], MIE = outA$MIE[i] + rnorm(1, mean = 0, sd = std),
                                             KE1 = outA$KE1[i] + rnorm(1, mean = 0, sd = std), 
                                             KE2 = outA$KE2[i] + rnorm(1, mean = 0, sd = std), REPLICATE = j))
  }
}
write.csv(fakedataA, paste0(output_folder, "/fakedata_qAOPA.csv"))

# Initial state and parameters for model C
inistateC = c(MIE = 0.5, KE1 = 0, KE2 = 0)
finishC <- seq(1, 100, by = 5)
parsC = c(tau = 0.05, eps = 0.2, d1 = 0.5, k_12 = 0.3, d2 = 0.4, k_2 = 0.2, h_2 = 0.3)
outC = as.data.frame(ode(inistateC, finishC, func = qAOPC, parms = parsC))

std <- 0.3  # standard deviation to add noise to the solution

fakedataC <- data.frame(time = numeric(), MIE = numeric(), KE1 = numeric(), KE2 = numeric(), REPLICATE = numeric())

for (i in 1:length(finishC)) {
  for (j in 1:nrep) {
    fakedataC <- rbind(fakedataC, data.frame(time = finishC[i], MIE = outC$MIE[i] + rnorm(1, mean = 0, sd = std),
                                             KE1 = outC$KE1[i] + rnorm(1, mean = 0, sd = std), 
                                             KE2 = outC$KE2[i] + rnorm(1, mean = 0, sd = std), REPLICATE = j))
  }
}
write.csv(fakedataC, paste0(output_folder, "/fakedata_qAOPC.csv"))

# Process the data for each model
data_stanA = fakedataA %>% 
  pivot_longer(!c(time, REPLICATE), names_to = "StateVar", values_to = "value") %>%
  group_by(time, StateVar) %>%
  summarise(mean = mean(value), sd = sd(value)) %>%
  ungroup()

data_stanB = fakedataB %>% 
  pivot_longer(!c(time, REPLICATE), names_to = "StateVar", values_to = "value") %>%
  group_by(time, StateVar) %>%
  summarise(mean = mean(value), sd = sd(value)) %>%
  ungroup()

data_stanC = fakedataC %>% 
  pivot_longer(!c(time, REPLICATE), names_to = "StateVar", values_to = "value") %>%
  group_by(time, StateVar) %>%
  summarise(mean = mean(value), sd = sd(value)) %>%
  ungroup()

# Define the input data for the cmdstan run
data_listA = list(
  N = length(unique(fakedataA$time)),
  t0 = 0,
  ts = unique(fakedataA$time), # does not have to include t0
  y_dose1 = as.matrix(data_stanA %>% select(-sd) %>%
                        group_by(time) %>%
                        pivot_wider(names_from = StateVar, values_from = mean) %>%
                        ungroup() %>% select(-time))[, c(3, 1, 2)],
  sigma_dose1 = as.matrix(data_stanA %>% select(-mean) %>%
                            group_by(time) %>%
                            pivot_wider(names_from = StateVar, values_from = sd) %>%
                            ungroup() %>% select(-time))[, c(3, 1, 2)],
  KE10 = 0,
  KE20 = 0
)

data_listB = list(
  N = length(unique(fakedataB$time)),
  t0 = 0,
  ts = unique(fakedataB$time), # does not have to include t0
  y_dose1 = as.matrix(data_stanB %>% select(-sd) %>%
                        group_by(time) %>%
                        pivot_wider(names_from = StateVar, values_from = mean) %>%
                        ungroup() %>% select(-time))[, c(3, 1, 2)],
  sigma_dose1 = as.matrix(data_stanB %>% select(-mean) %>%
                            group_by(time) %>%
                            pivot_wider(names_from = StateVar, values_from = sd) %>%
                            ungroup() %>% select(-time))[, c(3, 1, 2)],
  KE10 = 0,
  KE20 = 0
)

data_listC = list(
  N = length(unique(fakedataC$time)),
  t0 = 0,
  ts = unique(fakedataC$time), # does not have to include t0
  y_dose1 = as.matrix(data_stanC %>% select(-sd) %>%
                        group_by(time) %>%
                        pivot_wider(names_from = StateVar, values_from = mean) %>%
                        ungroup() %>% select(-time))[, c(3, 1, 2)],
  sigma_dose1 = as.matrix(data_stanC %>% select(-mean) %>%
                            group_by(time) %>%
                            pivot_wider(names_from = StateVar, values_from = sd) %>%
                            ungroup() %>% select(-time))[, c(3, 1, 2)],
  KE10 = 0,
  KE20 = 0
)
