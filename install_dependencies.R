#!/usr/bin/env Rscript
#
# GoldMiner R Package Installer
# This script automatically installs all required R packages for GoldMiner
#

cat("=====================================\n")
cat("GoldMiner R Package Installer\n")
cat("=====================================\n\n")

# Required packages for GoldMiner
packages <- c("igraph", "reshape2", "dplyr", "data.table",
              "optparse", "ggplot2", "Cairo")

# Check CRAN mirror
if (is.na(getOption("repos")[[1]])) {
  options(repos = "https://cloud.r-project.org/")
  cat("Setting CRAN repository to: https://cloud.r-project.org/\n\n")
}

cat("[*] Checking installed packages...\n")

# Check which packages are already installed
installed <- packages[packages %in% installed.packages()[,"Package"]]
missing <- packages[!(packages %in% installed.packages()[,"Package"])]

if (length(installed) > 0) {
  cat("✓ Already installed packages:", paste(installed, collapse=", "), "\n")
}

if (length(missing) > 0) {
  cat("✗ Missing packages:", paste(missing, collapse=", "), "\n\n")

  cat("[*] Installing missing packages...\n")
  cat("This may take several minutes.\n\n")

  # Install missing packages
  install.packages(missing, dependencies = TRUE)

  # Verify installation
  cat("\n[*] Verifying installation...\n")
  still_missing <- missing[!(missing %in% installed.packages()[,"Package"])]

  if (length(still_missing) == 0) {
    cat("✓ All packages installed successfully!\n")
  } else {
    cat("✗ Failed to install:", paste(still_missing, collapse=", "), "\n")
    cat("\nPlease try installing manually in R console:\n")
    cat("install.packages(c(", paste(shQuote(still_missing), collapse=", "), "))\n")
  }
} else {
  cat("✓ All packages are already installed!\n")
}

cat("\n=====================================\n")
cat("Installation complete!\n")
cat("You can now use GoldMiner.\n")
cat("=====================================\n")
