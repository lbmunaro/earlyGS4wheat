# Clean workspace
rm(list = objects())

# Packages ----
library(tidyverse)
library(asreml)

load('Data/ILYT_Blues-Gmatrix18.07.25.RData')

str(ILYT_blues)

ILYT_blues_filtered <- ILYT_blues |>
  filter(str_detect(Env, paste(c('Adv', 'Stp', 'Neo', 'Urb'), collapse = '|'))) |>
  droplevels()

str(ILYT_blues_filtered)

# Model parameters
k <- 5
model_name <- 'BLUP_YT2sNFA5.asr'

# BLUP
BLUP_YT2sNFA5.asr <- asreml(
  fixed = Blue ~ TraitEnv,
  random = ~ rr(TraitEnv, k):Gen + 
    diag(TraitEnv):Gen,
  weights = Wt,
  family = asr_gaussian(dispersion = 1),
  data = ILYT_blues_filtered,
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '8gb'
)

# Iterative updates for convergence
count <- 0

save(BLUP_YT2sNFA5.asr, ILYT_blues_filtered, file = paste0('Data/', model_name, '.RData'))

while (!isTRUE(BLUP_YT2sNFA5.asr$converge) && count < 100) {
  count <- count + 1
  BLUP_YT2sNFA5.asr <- update(BLUP_YT2sNFA5.asr)
  message('↻ Update ', count, ' (Convergence = ', BLUP_YT2sNFA5.asr$converge, ')')
  save(BLUP_YT2sNFA5.asr, ILYT_blues_filtered0, file = paste0('Data/', model_name, '.RData'))
}

# Report convergence status
if (isTRUE(BLUP_YT2sNFA5.asr$converge)) {
  message('✅ ', model_name, ' converged after ', count, ' update(s).')
} else {
  message('⚠️ ', model_name, ' did not converge after ', count, ' update(s).')
}