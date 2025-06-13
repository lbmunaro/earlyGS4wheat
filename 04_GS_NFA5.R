# Multi-Trait Multi-Environment ----
# This script fits a multi-trait, multi-environment genomic prediction model
# using ASReml for various phenotype subsets defined from IL_Pheno.
# It is designed for parallel execution in a SLURM array job.

# Clean workspace
rm(list = objects())  # Removes all objects from the environment to start fresh.

# Load required packages
library(tidyverse)    # Collection of R packages for data manipulation and visualization.
library(asreml)       # ASReml-R package for linear mixed model fitting.

# Handle SLURM array task ID passed as command-line argument
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])  # Convert the argument to an integer

# Load data ----
IL_Pheno <- readRDS('Data/IL_Pheno.rds')
IL.MT_Pheno <- readRDS('Data/IL.MT_Pheno.rds')
Ginv <- readRDS('Data/Ginv.rds')

# Create named list of filtered datasets for model training
df <- list(
  ILYT_2422 = IL_Pheno |> filter(Stage == 'YT') |> droplevels(),
  ILYT_22 = IL_Pheno |> filter(Stage == 'YT', Year == '2022') |> droplevels(),
  ILYT_23 = IL_Pheno |> filter(Stage == 'YT', Year == '2023') |> droplevels(),
  ILYT_2322 = IL_Pheno |> filter(Stage == 'YT', Year %in% c('2022', '2023')) |> droplevels(),
  IL_22 = IL_Pheno |> filter(Year == '2022') |> droplevels(),
  IL_23 = IL_Pheno |> filter(Year == '2023') |> droplevels(),
  IL_2322 = IL_Pheno |> filter(Year %in% c('2022', '2023'), Env != '22-SI') |> droplevels(),
  IL.MT_22 = IL.MT_Pheno |> filter(Year == '2022') |> droplevels(),
  IL.MT_23 = IL.MT_Pheno |> filter(Year == '2023') |> droplevels(),
  IL.MT_2322 = IL.MT_Pheno |> filter(Year %in% c('2022', '2023'), Env != '22-SI') |> droplevels()
)

# Apply preprocessing to each dataset in the list
# Add standardized phenotype and summary statistics
df <- map(df, function(d) {
  d %>%
    mutate(
      Pheno_z = scale(Pheno, center = TRUE, scale = TRUE)[,1],
      Pheno_mean = mean(Pheno, na.rm = TRUE),
      Pheno_sd = sd(Pheno, na.rm = TRUE)
    )
})

# Set model parameters ----
k <- 5
workspace <- '96gb'
df_names <- names(df)
model_name <- paste0(df_names[task_id], '_NFA', k, 'asr')
data_i <- df[[task_id]] |> droplevels()  # Subset data for this specific SLURM task

# Fit ASReml model ----
fit <- tryCatch({
  model <- asreml(
    fixed = Pheno_z ~ TraitEnv,
    random = as.formula(bquote(
      ~ rr(TraitEnv, .(k)):vm(Gkeep, Ginv) +
        diag(TraitEnv):vm(Gkeep, Ginv) +
        diag(TraitEnv):ide(Gkeep) +
        diag(TraitEnv):Block
    )),
    residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
    sparse = ~ TraitEnv:Gdrop,
    data = data_i,
    na.action = na.method(x = 'include'),
    maxit = 13,
    workspace = workspace
  )
  
  # Check convergence and update if necessary (up to 100 times)
  count <- 0
  while (!isTRUE(model$converge) && count < 100) {
    count <- count + 1
    model <- update(model)
    message('↻ Update ', count, ' (Convergence = ', model$converge, ')')
  }
  
  # Report model convergence status
  if (isTRUE(model$converge)) {
    message('✅ ', model_name, ' converged after ', count, ' update(s).')
  } else {
    message('⚠️ ', model_name, ' did not converge after ', count, ' update(s).')
  }
  
  # Save model object
  assign(model_name, model)
  save(list = model_name, file = paste0('Data/', model_name, '.RData'))
}, error = function(e) {
  # Handle model fitting errors
  warning('❌ Model failed for ', model_name, ' → ', e$message)
})
