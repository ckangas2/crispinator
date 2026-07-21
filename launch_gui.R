# launch_gui.R
# Purpose: Ask user for inputs -> Save Config -> Run Targets

library(tcltk)
library(targets)

# 1. Ask for Input File (MAGeCK Output)
# Filters for .txt or .csv files to be helpful
input_file <- tclvalue(tkgetOpenFile(
  title = "Select MAGeCK Gene Summary File",
  filetypes = "{{Text Files} {.txt .csv}} {{All Files} *}"
))

# Defensive: Stop if user hit "Cancel"
if (!nzchar(input_file)) stop("No input file selected. Operation cancelled.")

# 2. Ask for Output Folder
output_dir <- tclvalue(tkchooseDirectory(
  title = "Select Output Folder for Figures"
))

if (!nzchar(output_dir)) stop("No output directory selected. Operation cancelled.")

# 3. Save Configuration for _targets.R to read
# We save this as an RDS so the pipeline allows dynamic paths
config <- list(
  input_path = input_file,
  output_dir = output_dir
)
saveRDS(config, "_targets_config.rds")

message("Configuration saved.")
message("Input: ", input_file)
message("Output: ", output_dir)

# 4. Run the Pipeline
message("\n--- Starting Pipeline ---\n")
tar_make()

message("\n--- Pipeline Complete ---")