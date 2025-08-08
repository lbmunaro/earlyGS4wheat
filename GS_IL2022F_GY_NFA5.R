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

str(IL_blues)

# Vector with genotypes present in 24YT and SI22
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

# Dataset with different scenarios
# Filter dataset to include only the four locations of interest and SI
IL_blues_filtered <- IL_blues |>
  filter(str_detect(Env, paste(c('Adv', 'Stp', 'Neo', 'Urb', 'SI'), collapse = '|'))) |>
  droplevels()
str(IL_blues_filtered)

## Genotypes from 'IL2022F' cohort that are in SI-23 are masked
data_list <- list(
  ########## GY ##########
  
  # 1: GY 22-23YT
  GmaskIL2022F_GY2223 = IL_blues_filtered |>
    filter(Trait=='GY') |>
    filter(Env %in% c('22-Adv','22-Neo','22-Stp','22-Urb',
                      '23-Adv','23-Neo','23-Stp','23-Urb')) |>
    droplevels(),
  
  # 2: GY 22-24YTnTRN
  GmaskIL2022F_GY2224nTRN = IL_blues_filtered |>
    filter(Trait=='GY') |>
    filter(Env %in% c('22-Adv','22-Neo','22-Stp','22-Urb',
                      '23-Adv','23-Neo','23-Stp','23-Urb',
                      '24-Adv','24-Neo','24-Stp','24-Urb')) |>
    filter(!Gkeep %in% IL2022Fcohort) |>
    droplevels(),
  
  # 3: GY 22-24YTnTRN + IL2022F
  GmaskIL2022F_GY2224nTRN_IL2022F = IL_blues_filtered |>
    filter(Trait=='GY') |>
    filter(
      (Env %in% c('22-Adv','22-Neo','22-Stp','22-Urb',
                  '23-Adv','23-Neo','23-Stp','23-Urb',
                  '24-Adv','24-Neo','24-Stp','24-Urb') & !Gkeep %in% IL2022Fcohort) |
        (Env == '23-SI' & Gkeep %in% IL2022Fcohort)
    ) |>
    droplevels(),
  
  # 4: GY 22-24YTyTRN
  GmaskIL2022F_GY2224yTRN = IL_blues_filtered |>
    filter(Trait=='GY') |>
    filter(
      Env %in% c('22-Adv','22-Neo','22-Stp','22-Urb',
                 '23-Adv','23-Neo','23-Stp','23-Urb',
                 '24-Urb') |
        (Env %in% c('24-Adv','24-Neo','24-Stp') & !Gkeep %in% IL2022Fcohort)
    ) |>
    droplevels(),
  
  # 5: GY 22-24YTyTRN + IL2022F
  GmaskIL2022F_GY2224yTRN_IL2022F = IL_blues_filtered |>
    filter(Trait=='GY') |>
    filter(
      Env %in% c('22-Adv','22-Neo','22-Stp','22-Urb',
                 '23-Adv','23-Neo','23-Stp','23-Urb',
                 '24-Urb') |
        (Env %in% c('24-Adv','24-Neo','24-Stp') & !Gkeep %in% IL2022Fcohort) |
        (Env == '23-SI' & Gkeep %in% IL2022Fcohort)
    ) |>
    droplevels()
)

# Set model parameters
k <- 5
workspace <- '48gb'
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
      ~ rr(TraitEnv, .(k)):vm(Gkeep, Ginv.sparse) +
        diag(TraitEnv):vm(Gkeep, Ginv.sparse) +
        diag(TraitEnv):ide(Gkeep)
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
  save(list = c(model_name, dataset_name, 'Ginv.sparse'), file = model_file)
  message("✅ Initial model saved: ", model_file)
  
  # Iteratively update model until convergence or 100 updates
  while (!isTRUE(model$converge) && count < 100) {
    count <- count + 1
    model <- update(model)
    message('↻ Update ', count, ' (Convergence = ', model$converge, ')')
    assign(model_name, model)
    assign(dataset_name, data_i)
    save(list = c(model_name, dataset_name, 'Ginv.sparse'), file = model_file)
    message("💾 Model saved after update ", count)
  }
  
  # 🔧 Final save (post-loop)
  assign(model_name, model)
  assign(dataset_name, data_i)
  save(list = c(model_name, dataset_name, 'Ginv.sparse'), file = model_file)
  message("✔ Final model saved after convergence (or reaching iteration limit).")
  
}, error = function(e) {
  message("❌ Error during model fitting: ", conditionMessage(e))
  traceback(2)
})