#!/bin/bash
#SBATCH --mail-user=jzhan367@jh.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cpel_old_young
#SBATCH --chdir=.
#SBATCH --partition=shared
#SBATCH --mem=200G
#SBATCH --time=400:00:00
#SBATCH --output=cpel_old_young.out
#SBATCH --error=cpel_old_young.err
#SBATCH --cpus-per-task=8

# Get current date in yy_mm_dd_hh_mm_ss format
DATE=$(date +%y_%m_%d_%H_%M_%S)

# Create main output directory
MAIN_DIR="cpel_old_young_${DATE}"
mkdir -p "$MAIN_DIR"

echo "Created main directory: $MAIN_DIR"

# Loop through chromosomes 1-19
for chr in {1..19}; do
    # Create chromosome subdirectory
    CHR_DIR="${MAIN_DIR}/chr_${chr}"
    mkdir -p "$CHR_DIR"
    
    echo "Processing chromosome $chr"
    echo "Created subdirectory: $CHR_DIR"
    
    # Run Julia script with chromosome number and output directory
    julia cpel.jl "$CHR_DIR" "$chr"
    
    echo "Completed processing chromosome $chr"
done

echo "All chromosomes processed successfully!"

