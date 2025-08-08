# Clean workspace
rm(list = objects())

# Packages ----
library(tidyverse)
library(asreml)

load('Data/ILYT_Blues-Gmatrix18.07.25.RData')

str(ILYT_blues)

# Model parameters
k <- 5
model_name <- 'BLUP_2stpNFA5.asr'

# BLUP
BLUP_2stpNFA5.asr <- asreml(
  fixed = Blue ~ TraitEnv,
  random = ~ rr(TraitEnv, k):Gen + 
    diag(TraitEnv):Gen,
  weights = Wt,
  family = asr_gaussian(dispersion = 1),
  data = ILYT_blues,
  na.action = na.method(x = 'include'),
  maxit = 13,
  workspace = '8gb'
)

# Iterative updates for convergence
count <- 0

save(BLUP_2stpNFA5.asr, ILYT_blues, file = paste0('Data/', model_name, '.RData'))

while (!isTRUE(BLUP_2stpNFA5.asr$converge) && count < 100) {
  count <- count + 1
  BLUP_2stpNFA5.asr <- update(BLUP_2stpNFA5.asr)
  message('↻ Update ', count, ' (Convergence = ', BLUP_2stpNFA5.asr$converge, ')')
  save(BLUP_2stpNFA5.asr, ILYT_blues, file = paste0('Data/', model_name, '.RData'))
}

# Report convergence status
if (isTRUE(BLUP_2stpNFA5.asr$converge)) {
  message('✅ ', model_name, ' converged after ', count, ' update(s).')
} else {
  message('⚠️ ', model_name, ' did not converge after ', count, ' update(s).')
}