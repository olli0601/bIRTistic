#!/usr/bin/env Rscript
# Run all tests in the test directory

library(testthat)

# Set working directory to project root
if (grepl("test$", getwd())) {
  setwd("..")
} else if (!file.exists("test")) {
  stop("Cannot find test directory")
}

cat("Running all tests in test/ directory\n")
cat("=====================================\n\n")

# Run all test files
test_results <- test_dir("test", reporter = "progress")

cat("\n")
cat("=====================================\n")
cat("All tests completed\n")
