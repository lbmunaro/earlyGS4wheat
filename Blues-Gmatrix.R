# Clean workspace
rm(list = objects())

# Packages ----
library(tidyverse)
library(asreml)
library(future.apply)

load('Data/ILYT_Pheno-Gmatrix07.07.25.RData')

trait_map <- c(
  "Grain.yield...bu.ac" = 'GY',
  "Grain.test.weight...lbs.bu" = 'TW',
  "Heading.time...Julian.date" = 'HD',
  "Plant.height.inches" = 'HT',
  "Maturity.time" = 'MAT'
)

# Blues from YT - 2022 to 2025
ILYT_blues <- ILYT_Pheno |>
  filter(trait %in% names(trait_map)) |>
  mutate(Trait = recode(trait, !!!trait_map)) |>
  filter(conv==TRUE) |>
  mutate(
    study = str_replace(study,'YT_Stj_22','YT_Addie_22'),
    study = as.factor(gsub('^YT_', '', study)),
    study = as.factor(gsub('^Addie_', 'Adv_', study)),
    study = str_replace(study, '(\\w+)_(\\d+)', '\\2-\\1')
  ) |>
  rename(Env=study,
         Blue=Pheno,
         Wt = weight) |>
  mutate(TraitEnv=as.factor(paste(Trait,Env,sep = '-')),
         Gen = factor(if_else(str_detect(Gen, 'FILL'), NA, Gen)),
         Gdrop = factor(if_else(str_detect(Gdrop, 'FILL'), NA, Gdrop))
         ) |>
  select(TraitEnv, Trait, Env, Gen, Gkeep, Gdrop, rel, Blue, Wt) |>
  mutate_if(is.character, ~as.factor(.x)) |>
  droplevels()

# Phenotypic data from SI - 2022 and 2023
ILSI_Pheno <- readRDS('Data/IL.MT_Pheno.rds') |>
  filter(str_detect(Env,'-SI')) |>
  mutate(Pheno = ifelse(Trait=='GY',Pheno/67.251,Pheno),
         Pheno = ifelse(Trait=='TW',Pheno/1000*2.2046*35.2391,Pheno)) |>
  mutate(
    Gkeep = factor(if_else(Gen %in% attr(Ginv.sparse, "rowNames"), Gen, NA)),  # Genotypes with marker data.
    Gdrop = factor(if_else(Gen %in% attr(Ginv.sparse, "rowNames"), NA, Gen))  # Genotypes without marker data.
  ) |>
  filter(TraitEnv!='HT-24-SI1') |>
  droplevels()

# Calculate wt and rel for SI

## BLUP model ----
# Function to calculate rel and Wt, this model will run for each TraitEnv in ILSI_Pheno
mod.blups <- function(dat){
  mod <- asreml(fixed = Pheno~1,
                random = ~vm(Gkeep, Ginv.sparse),
                residual = ~ar1(Col):ar1(Row),
                data = dat,
                na.action = na.method(x='include'),
                maxit = 13,
                workspace = '1gb')
  # Update model
  mod <- update(mod)
  # Extract BLUPs and calculate reliability
  blups<- predict(mod, classify='Gkeep', ignore=c('(Intercept)'),
                  pworkspace='2gb')$pvals
  pev<- blups[,'std.error']^2 # Prediction error variance
  Vg<- summary(mod)$varcomp['vm(Gkeep, Ginv.sparse)','component'] # Genetic variance estimate
  rel<- 1-(pev/Vg) # Reliability for each genotype
#  rel<- mean(rel) # Mean reliability
  # Return mean reliability
  res<- data.frame(Gkeep = blups[['Gkeep']],
                   rel=rel,
                   Wt=(1/pev))
  return(res)
}

# Plan parallel execution
plan(multisession, workers = availableCores() - 1)

# Prepare data: split by TraitEnv
traitenv_list <- split(ILSI_Pheno, ILSI_Pheno$TraitEnv)

# Define parallel-safe wrapper
run_blup_model <- function(dat) {
  result <- tryCatch(mod.blups(dat), error = function(e) {
    message("Failed for ", unique(dat$TraitEnv), ": ", e$message)
    return(NULL)
  })
  
  if (is.null(result)) return(NULL)
  
  dat_out <- dat |>
    select(TraitEnv, Trait, Env, Gen, Gkeep, Gdrop, Pheno) |>
    rename(Blue = Pheno) |>
    left_join(result, by = "Gkeep")  # merges rel and Wt
  
  return(dat_out)
}

# Run models in parallel
ILSI_blues_list <- future_lapply(traitenv_list, run_blup_model)

# Combine results into one data.frame
ILSI_blues <- bind_rows(ILSI_blues_list) |>
  mutate_if(is.character, as.factor) |>
  droplevels() |>
  group_by(TraitEnv) |>
  mutate(rel=mean(rel, na.rm=T)) |>
  ungroup()

str(ILSI_blues)

IL_blues <- bind_rows(ILYT_blues, ILSI_blues)

save(Ginv.sparse, ILYT_blues, ILSI_blues, IL_blues,
     file = 'Data/ILYT_Blues-Gmatrix07.22.25.RData')

# 
load('Data/ILYT_Blues-Gmatrix07.22.25.RData')

ILSI_blues |> 
  filter(Env=='24-SI1') |>
  filter(!is.na(Gkeep)) |>
  group_by(TraitEnv) |> reframe(rel=unique(rel))
