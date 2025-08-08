# Clean workspace
rm(list = objects())

# Packages ----
library(tidyverse)
library(asreml)

# Handle SLURM array task ID passed as command-line argument
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])  # Convert the argument to an integer

# Load data
load('Data/ILYT_Blues-Gmatrix07.21.25.RData')

# Vector with genotypes present in 24YT and SI22
SI23cohort <- ILSI_blues |>
  filter(str_detect(Env,'23-')) |>
  filter(str_starts(Gkeep, 'IL2022-')) |>
  droplevels() |>
  group_by(Gkeep) |>
  reframe() |>
  pull() |> as.vector()

SI23YTcohort <- ILYT_blues |>
  filter(str_detect(Env,'24-')) |>
  filter(str_starts(Gkeep, 'IL2022-')) |>
  droplevels() |>
  group_by(Gkeep) |>
  reframe() |>
  pull() |> as.vector()

# Vector with genotypes present in 25YT and SI23
IL2022Fcohort <- ILSI_blues |>
  filter(str_detect(Env,'23-')) |>
  filter(str_starts(Gkeep, 'IL2022F')) |>
  droplevels() |>
  group_by(Gkeep) |>
  reframe() |>
  pull() |> as.vector()

IL2022FYTcohort <- ILYT_blues |>
  filter(str_detect(Env,'25-')) |>
  filter(str_starts(Gkeep, 'IL2022F')) |>
  droplevels() |>
  group_by(Gkeep) |>
  reframe() |>
  pull() |> as.vector()

# Vector with genotypes present in 24YT and SI24
SI24cohort <- ILSI_blues |>
  filter(str_detect(Env,'24-')) |>
  filter(str_starts(Gkeep, 'IL23-')) |>
  droplevels() |>
  group_by(Gkeep) |>
  reframe() |>
  pull() |> as.vector()

SI24YTcohort <- ILYT_blues |>
  filter(str_detect(Env,'25-')) |>
  filter(str_starts(Gkeep, 'IL23-')) |>
  droplevels() |>
  group_by(Gkeep) |>
  reframe() |>
  pull() |> as.vector()

# Dataset with different scenarios
# Filter dataset to include only the four locations of interest and SI
IL_blues_filtered <- IL_blues |>
  filter(str_detect(Env, paste(c('Adv', 'Stp', 'Neo', 'Urb', 'SI'), collapse = '|'))) |>
  droplevels()
str(IL_blues_filtered)

## Genotypes from 'IL23-' cohort that are in SI-24 are masked
data_list <- list(
  
  #1: GY 22-23YTyTRN + SI23 - SI23
  BlupAH3SI23_GY2223yTRN = IL_blues_filtered |>
    filter(Trait=='GY') |>
    filter(Env%in%c('22-Adv','22-Neo','22-Stp','22-Urb',
                    '23-Adv','23-Neo','23-Stp','23-Urb')) |>
    droplevels(),
  
  #2: MT 22-23YTyTRN - SI23
  BlupAH3SI23_MT2223yTRN = IL_blues_filtered |>
    filter(Env%in%c('22-Adv','22-Neo','22-Stp','22-Urb',
                    '23-Adv','23-Neo','23-Stp','23-Urb')) |>
    droplevels(),
  
  #3: GY 22-24YTyTRN
  BlupAH3IL2022F_GY2224yTRN = IL_blues_filtered |>
    filter(Trait=='GY') |>
    filter(
      Env %in% c('22-Adv','22-Neo','22-Stp','22-Urb',
                 '23-Adv','23-Neo','23-Stp','23-Urb',
                 '24-Urb') |
        (Env %in% c('24-Adv','24-Neo','24-Stp') & !Gkeep %in% IL2022Fcohort)
    ) |>
    droplevels(),
  
  #4: MT 22-24YTyTRN
  BlupAH3IL2022F_MT2224yTRN = IL_blues_filtered |>
    filter(
      Env %in% c('22-Adv','22-Neo','22-Stp','22-Urb',
                 '23-Adv','23-Neo','23-Stp','23-Urb',
                 '24-Urb') |
        (Env %in% c('24-Adv','24-Neo','24-Stp') & !Gkeep %in% IL2022Fcohort)
    ) |>
    droplevels(),
  
  #5: GY 22-25YTyTRN
  BlupAH3SI24_GY2225yTRN = IL_blues_filtered |>
    filter(Trait=='GY') |>
    filter(Env %in% c('22-Adv','22-Neo','22-Stp','22-Urb',
                      '23-Adv','23-Neo','23-Stp','23-Urb',
                      '24-Adv','24-Neo','24-Stp','24-Urb',
                      '25-Adv','25-Neo','25-Stp','25-Urb')) |>
    droplevels(),
  
  #6: MT 22-25YTyTRN
  BlupAH3SI24_MT2225yTRN = IL_blues_filtered |>
    filter(Env %in% c('22-Adv','22-Neo','22-Stp','22-Urb',
                      '23-Adv','23-Neo','23-Stp','23-Urb',
                      '24-Adv','24-Neo','24-Stp','24-Urb',
                      '25-Adv','25-Neo','25-Stp','25-Urb')) |>
    droplevels()
)

# Set model parameters
k <- 5
workspace <- '8gb'
df_names <- names(data_list)
dataset_name <- df_names[task_id]
data_i <- data_list[[task_id]] |> droplevels()
model_name <- paste0(dataset_name, '_NFA', k)
model_file <- paste0('Data/', model_name, '.RData')

# Fit ASReml model ----
fit <- tryCatch({
  model <- asreml(
    fixed = Blue ~ TraitEnv,
    random = as.formula(bquote(
      ~ rr(TraitEnv, .(k)):Gkeep +
        diag(TraitEnv):Gkeep
    )),
    weights = Wt,
    family = asr_gaussian(dispersion = 1),
    data = data_i,
    na.action = na.method(x = 'include'),
    maxit = 13,
    workspace = workspace
  )
  
  # Save initial model state
  count <- 0
  assign(model_name, model)
  assign(dataset_name, data_i)
  save(list = c(model_name, dataset_name), file = model_file)
  message("✅ Initial model saved: ", model_file)
  
  # Iteratively update model until convergence or 100 updates
  while (!isTRUE(model$converge) && count < 100) {
    count <- count + 1
    model <- update(model)
    message('↻ Update ', count, ' (Convergence = ', model$converge, ')')
    assign(model_name, model)
    assign(dataset_name, data_i)
    save(list = c(model_name, dataset_name), file = model_file)
    message("💾 Model saved after update ", count)
  }
  
  # 🔧 Final save (post-loop)
  assign(model_name, model)
  assign(dataset_name, data_i)
  save(list = c(model_name, dataset_name), file = model_file)
  message("✔ Final model saved after convergence (or reaching iteration limit).")
  
}, error = function(e) {
  message("❌ Error during model fitting: ", conditionMessage(e))
  traceback(2)
})