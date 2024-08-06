if (!require("pacman", quietly = T)) {
  install.packages("pacman")
}
pacman::p_load(tidyverse, cmdstanr, posterior, bayesplot, deSolve, conflicted,
               loo, gridExtra, latex2exp,viridis,scales)

conflict_prefer("sd", "stats")
conflict_prefer("filter","dplyr")
source("~/Phd material/qAOP_framework/model_comparison/qAOP_models.R")
source("~/Phd material/qAOP_framework/model_comparison/fakedataABC.R")

model = readline("choose model to run: ")
data = readline("chose data to use: ")
data_list = get(paste0("data_list",data))

options(mc.cores = parallel::detectCores())  # Use all available cores
#, stanc.allow_optimizations = TRUE, stanc.auto_format = TRUE
chains <- 4 # Number of chains
iter <- 10000  # Number of iterations
warmup <- 6000 # Number of warmup iterations
thin <- 1  # Thinning parameter

compiled_model <- cmdstan_model(paste0("/home/ditilliof/Stan_models/qAOP",model,".stan")) #modify according to you file path

# This is for cmdstanr
fit <- compiled_model$sample(data = data_list, 
                             parallel_chains = getOption("mc.cores", 4), 
                             chains=chains,
                             iter_warmup = warmup, 
                             iter_sampling = iter-warmup,
                             #show_messages = FALSE,
                             refresh = 1
)
parsA = c("tau", "eps", "d1", "k_12", "d2", "k_2", "S0")
parsB = c("tau", "eps", "d1", "k_12", "d2", "S0")
parsC = c("tau", "eps", "d1", "k_12", "d2", "k_2", "h_2", "S0")
pars = get(paste0("pars",model))
# Extract draws and convert to data frame with iteration and chain
draws_df <- fit$draws(format = "df", variables = pars) %>%
  mutate(.iteration = rep(1:(nrow(.) / max(.chain)), times = max(.chain)),
         .chain = rep(1:max(.chain), each = nrow(.) / max(.chain)))

# Custom labels using expressions for ggplot2
custom_labelsA <- c(tau = "tau", eps = "epsilon", d1 = "d[1]", k_12 = "k[12]", d2 = "d[2]",
                    k_2 = "k[2]", S0 = "S[0]")
custom_labelsB <- c(tau = "tau", eps = "epsilon", d1 = "d[1]", k_12 = "k[12]", d2 = "d[2]",
  S0 = "S[0]")
custom_labelsC <- c(tau = "tau", eps = "epsilon", d1 = "d[1]", k_12 = "k[12]", d2 = "d[2]",
                    k_2 = "k[2]", h_2 = "h[2]", S0 = "S[0]")
custom_labels = get(paste0("custom_labels", model))

# Remove unnecessary columns
draws_df <- draws_df %>% select(-.draw)

# Base theme with larger axis title, facet text size, legend text size, and tick labels
base_theme <- theme_minimal() + 
  theme(
    axis.title = element_text(size = 18),  # Larger axis titles
    strip.text = element_text(size = 18),  # Larger facet labels
    legend.title = element_text(size = 18), # Larger legend title
    legend.text = element_text(size = 18),  # Larger legend text
    axis.text = element_text(size = 14)     # Larger tick labels
  )

# Function to round axis labels to 3 decimal places
round3 <- label_number(accuracy = 0.001)

# Create trace plot
fit_trace <- draws_df %>%
  pivot_longer(cols = -c(.iteration, .chain), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = .iteration, y = value, color = as.factor(.chain))) +
  geom_line(alpha = 0.5) +
  facet_wrap(~variable, scales = "free", labeller = labeller(variable = as_labeller(custom_labels, label_parsed))) +
  scale_color_viridis_d() +
  labs(x = "Iteration", y = "Parameter value", color = "Chain") +
  scale_x_continuous(breaks = c(0, 2000, 4000)) +  # Adjust x-axis breaks
  base_theme

# Create density plot
fit_density <- draws_df %>%
  pivot_longer(cols = -c(.iteration, .chain), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, fill = as.factor(.chain))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~variable, scales = "free", labeller = labeller(variable = as_labeller(custom_labels, label_parsed))) +
  scale_fill_viridis_d() +
  labs(x = "Parameter value", y = "Density", fill = "Chain") +
  scale_x_continuous(labels = round3, breaks = function(x) pretty(x, n = 3)) +  # Adjust x-axis breaks and round labels
  base_theme

# Create density plot with vertical lines
# fit_density_lines <- draws_df %>%
#   pivot_longer(cols = -c(.iteration, .chain), names_to = "variable", values_to = "value") %>%
#   ggplot(aes(x = value, fill = as.factor(.chain))) +
#   geom_density(alpha = 0.5) +
#   facet_wrap(~variable, scales = "free", labeller = labeller(variable = as_labeller(custom_labels, label_parsed))) +
#   scale_fill_viridis_d() +
#   labs(x = "Parameter value", y = "Density", fill = "Chain") +
#   geom_vline(data = data.frame(variable = c("tau", "eps", "d1", "k_12", "d2", "S0"), xintercept = c(0.05, 0.2, 0.5, 0.3, 0.4, 0.5)), 
#              aes(xintercept = xintercept), color = "red", linetype = "dashed") +
#   scale_x_continuous(labels = round3, breaks = function(x) pretty(x, n = 3)) +  # Adjust x-axis breaks and round labels
#   base_theme
#pars = c(tau=0.05, eps = 0.2, d1 = 0.5, k_12 = 0.3, d2 = 0.4)
# Save the plots as TIFF files
tiff(paste0(output_folder, "/", "Trace-plot_qAOP", model,"_datafrom",data, ".png"), units = "in", width = 9, height = 5.3, res = 700, compression = 'lzw')
print(fit_trace)
dev.off()

tiff(paste0(output_folder, "/", "Density-plot_qAOP", model,"_datafrom",data,".png"), units = "in", width = 9, height = 5.3, res = 700, compression = 'lzw')
print(fit_density)
dev.off()

# tiff(paste0(output_folder, "/", "Density-plot_qAOPB_lines", ".png"), units = "in", width = 9, height = 5.3, res = 700, compression = 'lzw')
# print(fit_density_lines)
# dev.off()

# Display the plots together
grid.arrange(fit_trace, fit_density, ncol = 1)

# Save the draws and summaries
draws <- fit$draws(format = "draws_matrix", variables = pars)
fit_summary <- fit$summary()

loglikelihood <- fit$draws(format = "draws_matrix", variables = 'logLikelihood')
draws_complete <- fit$draws(format = "draws_matrix")

sumPar <- data.frame(fit_summary) %>% column_to_rownames("variable")

parSets <- draws

x = fit$loo(variables = "logLikelihood", cores = 3)
estimates = data.frame(x$estimates)
write.csv(estimates, paste0(output_folder, "/loo_output_qAOP", model,"_datafrom",data,".csv"))
saveRDS(fit, file = paste0(output_folder, "/", "fitqAOP", model,"_datafrom",data,".rds"))
write_csv(data.frame(fit_summary), paste0(output_folder, "/", "parameter_summary_qAOP", model,"_datafrom",data,".csv"))
write_csv(data.frame(draws_complete), paste0(output_folder, "/", "draws_qAOP", model,"_datafrom",data,".csv"))

model = readline("choose second model to run: ")

options(mc.cores = parallel::detectCores())  # Use all available cores
#, stanc.allow_optimizations = TRUE, stanc.auto_format = TRUE
chains <- 4 # Number of chains
iter <- 10000  # Number of iterations
warmup <- 6000 # Number of warmup iterations
thin <- 1  # Thinning parameter

compiled_model <- cmdstan_model(paste0("/home/ditilliof/Stan_models/qAOP",model,".stan")) #modify according to you file path

# This is for cmdstanr
fit <- compiled_model$sample(data = data_list, 
                             parallel_chains = getOption("mc.cores", 4), 
                             chains=chains,
                             iter_warmup = warmup, 
                             iter_sampling = iter-warmup,
                             #show_messages = FALSE,
                             refresh = 1,
                             #seed = 1234
)
parsA = c("tau", "eps", "d1", "k_12", "d2", "k_2", "S0")
parsB = c("tau", "eps", "d1", "k_12", "d2", "S0")
parsC = c("tau", "eps", "d1", "k_12", "d2", "k_2", "h_2", "S0")
pars = get(paste0("pars",model))
# Extract draws and convert to data frame with iteration and chain
draws_df <- fit$draws(format = "df", variables = pars) %>%
  mutate(.iteration = rep(1:(nrow(.) / max(.chain)), times = max(.chain)),
         .chain = rep(1:max(.chain), each = nrow(.) / max(.chain)))

# Custom labels using expressions for ggplot2
custom_labelsA <- c(tau = "tau", eps = "epsilon", d1 = "d[1]", k_12 = "k[12]", d2 = "d[2]",
                    k_2 = "k[2]", S0 = "S[0]")
custom_labelsB <- c(tau = "tau", eps = "epsilon", d1 = "d[1]", k_12 = "k[12]", d2 = "d[2]",
                    S0 = "S[0]")
custom_labelsC <- c(tau = "tau", eps = "epsilon", d1 = "d[1]", k_12 = "k[12]", d2 = "d[2]",
                    k_2 = "k[2]", h_2 = "h[2]", S0 = "S[0]")
custom_labels = get(paste0("custom_labels", model))

# Remove unnecessary columns
draws_df <- draws_df %>% select(-.draw)

# Base theme with larger axis title, facet text size, legend text size, and tick labels
base_theme <- theme_minimal() + 
  theme(
    axis.title = element_text(size = 18),  # Larger axis titles
    strip.text = element_text(size = 18),  # Larger facet labels
    legend.title = element_text(size = 18), # Larger legend title
    legend.text = element_text(size = 18),  # Larger legend text
    axis.text = element_text(size = 14)     # Larger tick labels
  )

# Function to round axis labels to 3 decimal places
round3 <- label_number(accuracy = 0.001)

# Create trace plot
fit_trace <- draws_df %>%
  pivot_longer(cols = -c(.iteration, .chain), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = .iteration, y = value, color = as.factor(.chain))) +
  geom_line(alpha = 0.5) +
  facet_wrap(~variable, scales = "free", labeller = labeller(variable = as_labeller(custom_labels, label_parsed))) +
  scale_color_viridis_d() +
  labs(x = "Iteration", y = "Parameter value", color = "Chain") +
  scale_x_continuous(breaks = c(0, 2000, 4000)) +  # Adjust x-axis breaks
  base_theme

# Create density plot
fit_density <- draws_df %>%
  pivot_longer(cols = -c(.iteration, .chain), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value, fill = as.factor(.chain))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~variable, scales = "free", labeller = labeller(variable = as_labeller(custom_labels, label_parsed))) +
  scale_fill_viridis_d() +
  labs(x = "Parameter value", y = "Density", fill = "Chain") +
  scale_x_continuous(labels = round3, breaks = function(x) pretty(x, n = 3)) +  # Adjust x-axis breaks and round labels
  base_theme

# Create density plot with vertical lines
# fit_density_lines <- draws_df %>%
#   pivot_longer(cols = -c(.iteration, .chain), names_to = "variable", values_to = "value") %>%
#   ggplot(aes(x = value, fill = as.factor(.chain))) +
#   geom_density(alpha = 0.5) +
#   facet_wrap(~variable, scales = "free", labeller = labeller(variable = as_labeller(custom_labels, label_parsed))) +
#   scale_fill_viridis_d() +
#   labs(x = "Parameter value", y = "Density", fill = "Chain") +
#   geom_vline(data = data.frame(variable = c("tau", "eps", "d1", "k_12", "d2", "S0"), xintercept = c(0.05, 0.2, 0.5, 0.3, 0.4, 0.5)), 
#              aes(xintercept = xintercept), color = "red", linetype = "dashed") +
#   scale_x_continuous(labels = round3, breaks = function(x) pretty(x, n = 3)) +  # Adjust x-axis breaks and round labels
#   base_theme
#pars = c(tau=0.05, eps = 0.2, d1 = 0.5, k_12 = 0.3, d2 = 0.4)
# Save the plots as TIFF files
tiff(paste0(output_folder, "/", "Trace-plot_qAOP", model,"_datafrom",data, ".png"), units = "in", width = 9, height = 5.3, res = 700, compression = 'lzw')
print(fit_trace)
dev.off()

tiff(paste0(output_folder, "/", "Density-plot_qAOP", model,"_datafrom",data,".png"), units = "in", width = 9, height = 5.3, res = 700, compression = 'lzw')
print(fit_density)
dev.off()

# tiff(paste0(output_folder, "/", "Density-plot_qAOPB_lines", ".png"), units = "in", width = 9, height = 5.3, res = 700, compression = 'lzw')
# print(fit_density_lines)
# dev.off()

# Display the plots together
grid.arrange(fit_trace, fit_density, ncol = 1)

# Save the draws and summaries
draws <- fit$draws(format = "draws_matrix", variables = pars)
fit_summary <- fit$summary()

loglikelihood <- fit$draws(format = "draws_matrix", variables = 'logLikelihood')
draws_complete <- fit$draws(format = "draws_matrix")

sumPar <- data.frame(fit_summary) %>% column_to_rownames("variable")

parSets <- draws
x = fit$loo(variables = "logLikelihood", cores = 3)
estimates = data.frame(x$estimates)
write.csv(estimates, paste0(output_folder, "/loo_output_qAOP", model,"_datafrom",data,".csv"))

saveRDS(fit, file = paste0(output_folder, "/", "fitqAOP", model,"_datafrom",data,".rds"))
write_csv(data.frame(fit_summary), paste0(output_folder, "/", "parameter_summary_qAOP", model,"_datafrom",data,".csv"))
write_csv(data.frame(draws_complete), paste0(output_folder, "/", "draws_qAOP", model,"_datafrom",data,".csv"))

