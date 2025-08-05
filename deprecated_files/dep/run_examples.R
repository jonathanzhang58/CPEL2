# Example usage of the consolidated partitioning script
# This script demonstrates how to run the partitioning.R script with different landmark types

# Example 1: Run with promoters landmark (equivalent to 250bp.R)
cat("Running with promoters landmark...\n")
system("Rscript partitioning.R 250 promoters")

# Example 2: Run with shores landmark (equivalent to 250bpshores.R)
cat("Running with shores landmark...\n")
system("Rscript partitioning.R 250 shores")

# Example 3: Run with islands landmark (equivalent to 250bpislands.R)
cat("Running with islands landmark...\n")
system("Rscript partitioning.R 250 islands")

# Example 4: Run with different tile width
cat("Running with 500bp tiles and shores landmark...\n")
system("Rscript partitioning.R 500 shores")

cat("All examples completed!\n") 