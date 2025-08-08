# Clean workspace
rm(list = objects())

source('Functions_earlyGS.R')

load('Data/BLUP_YT2sNFA5.asr.RData')

blup <- blup_asreml(BLUP_YT2sNFA5.asr,5,'TraitEnv','Gen') |>
  separate(TE_fct, c('Trait','Year','Loc'), sep='-', remove = F) |>
  filter(Loc%in%c('Adv','Neo','Stp','Urb')) |>
  filter(Year!='22'&Year!='23')

str(blup)

IL2022FYTcohort <- ILYT_blues_filtered |>
  filter(str_detect(Env,'25-')) |>
  filter(str_starts(Gkeep, 'IL2022F')) |>
  droplevels() |>
  group_by(Gkeep) |>
  reframe() |>
  pull() |> as.vector()

# Initialize lists to collect GTEacc and Gacc results
GTEacc_list <- list()
Gacc_list <- list()
Conv_list <- list()

# Define a named vector of models and file paths
model_info <- list(
  GmaskIL2022F_GY2223_NFA5 = 'Data/GmaskIL2022F_GY2223_NFA5.RData',
  GmaskIL2022F_GY2224nTRN_NFA5 = 'Data/GmaskIL2022F_GY2224nTRN_NFA5.RData',
  GmaskIL2022F_GY2224nTRN_IL2022F_NFA5 = 'Data/GmaskIL2022F_GY2224nTRN_IL2022F_NFA5.RData',
  GmaskIL2022F_GY2224yTRN_NFA5 = 'Data/GmaskIL2022F_GY2224yTRN_NFA5.RData',
  GmaskIL2022F_GY2224yTRN_IL2022F_NFA5 = 'Data/GmaskIL2022F_GY2224yTRN_IL2022F_NFA5.RData',
  
  GmaskIL2022F_MT2223_NFA5 = 'Data/GmaskIL2022F_MT2223_NFA5.RData',
  GmaskIL2022F_MT2224nTRN_NFA5 = 'Data/GmaskIL2022F_MT2224nTRN_NFA5.RData',
  GmaskIL2022F_MT2224nTRN_IL2022F_NFA5 = 'Data/GmaskIL2022F_MT2224nTRN_IL2022F_NFA5.RData',
  GmaskIL2022F_MT2224yTRN_NFA5 = 'Data/GmaskIL2022F_MT2224yTRN_NFA5.RData',
  GmaskIL2022F_MT2224yTRN_IL2022F_NFA5 = 'Data/GmaskIL2022F_MT2224yTRN_IL2022F_NFA5.RData'
)

# Loop over each model
for (model_name in names(model_info)) {
  
  # Load model
  load(model_info[[model_name]])
  
  # Convergence
  Conv_temp <- get(model_name)$converge
  # Extract gebv
  gebv <- gebvs_asreml(get(model_name), 5, 'TraitEnv') |>
    rename(gebv = blup) |>
    separate(TE_fct, c('Trait','Year','Loc'), sep = '-', remove = FALSE) |>
    filter(Loc %in% c('Adv', 'Neo', 'Stp', 'Urb'))
  
  # GTEacc
  # GTEacc_tmp <- gebv |>
  #   left_join(blup) |>
  #   filter(G_fct %in% IL2022FYTcohort) |>
  #   group_by(Trait) |>
  #   summarise(GTEcor = cor(gebv, blup, use = "complete.obs"), .groups = "drop") |>
  #   mutate(Model = model_name)
  
  # Gacc
  Gacc_tmp <- gebv |>
    filter(G_fct %in% IL2022FYTcohort) |>
    group_by(Trait, G_fct) |>
    summarise(G_gebv = mean(gebv), .groups = "drop") |>
    left_join(
      blup |>
        filter(G_fct %in% IL2022FYTcohort) |>
        group_by(Trait, G_fct) |>
        summarise(G_blup = mean(blup), .groups = "drop"),
      by = c("Trait", "G_fct")
    ) |>
    group_by(Trait) |>
    summarise(Gcor = cor(G_gebv, G_blup, use = "complete.obs"), .groups = "drop") |>
    mutate(Model = model_name)
  
  # Append results to lists
  #GTEacc_list[[model_name]] <- GTEacc_tmp
  Gacc_list[[model_name]] <- Gacc_tmp
  Conv_list[[model_name]] <- Conv_temp
}

# Combine all results into summary tables
Conv_summary <- bind_cols(Conv_list) |> pivot_longer(cols = everything(),names_to = 'model', values_to = 'Conv')
#GTEacc_summary <- bind_rows(GTEacc_list) |> select(Model, Trait, GTEcor) |> arrange(Trait)
Gacc_summary   <- bind_rows(Gacc_list)   |> select(Model, Trait, Gcor) |> arrange(Trait)

# Show results
print(Conv_summary)
#print(GTEacc_summary)
print(Gacc_summary)
