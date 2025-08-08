rm(list = objects())

library(asreml)

# Command line argument for SLURM array task ID
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])


# List of model base names (without .RData extension)
model_files <- c(
  "GmaskSI22_MT2224_NFA5",
  "GmaskSI24_MT22.24_NFA5",
  "GmaskSI24_MT2225nTRN_NFA5",
  "GmaskSI24_MT2225nTRN_SI24_NFA5",
  "GmaskSI24_MT2225yTRN_NFA5",
  "GmaskSI24_MT2225yTRN_SI24_NFA5"
)

model_name <- model_files[task_id]
model_file <- paste0("Data/", model_name, ".RData")

# Load model and related objects
load(model_file)

# Recover objects from environment
model <- get(model_name)

# Infer dataset name by searching for any object that starts with model_name minus "_NFA5"
base_prefix <- sub("_NFA[0-9]+$", "", model_name)
dataset_name <- grep(base_prefix, ls(), value = TRUE)[1]
data_i <- get(dataset_name)
k <- 5
workspace <- '64gb'

# Iterative update loop
count <- 0
while (!isTRUE(model$converge) && count < 100) {
  count <- count + 1
  model <- update(model)
  message('↻ Update ', count, ' (Convergence = ', model$converge, ')')
  
  assign(model_name, model)
  assign(dataset_name, data_i)
  save(list = c(model_name, dataset_name, "Ginv.sparse"), file = model_file)
  message("💾 Model saved after update ", count)
}

# Final save (redundant but safe)
assign(model_name, model)
assign(dataset_name, data_i)
save(list = c(model_name, dataset_name, "Ginv.sparse"), file = model_file)
message("✔ Final model saved after convergence (or reaching iteration limit).")