#!/bin/bash
#SBATCH --mail-user=jzhan367@jh.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cpel_old_young
#SBATCH --chdir=.
#SBATCH --partition=shared,cegs2,cee
#SBATCH --mem-per-cpu=15G
#SBATCH --time=40:00:00
#SBATCH --output=cpel_old_young.out
#SBATCH --error=cpel_old_young.err
#SBATCH --cpus-per-task=8


module load conda_R


# Get current date in yy_mm_dd_hh_mm_ss format
DATE=$(date +%y_%m_%d_%H_%M_%S)

# Create main output directory
MAIN_DIR="cpel_old_young_${DATE}"
mkdir -p "$MAIN_DIR"

echo "Created main directory: $MAIN_DIR"

# Store the absolute path before changing directories
ABSOLUTE_MAIN_DIR="$(pwd)/$MAIN_DIR"

cd $MAIN_DIR


Rscript /dcs07/afeinber/data/personal/jzhan/CPEL2/250bp.R

# Run ROIcut.R with the absolute main directory as argument
echo "Running ROIcut.R with directory: $ABSOLUTE_MAIN_DIR"
Rscript /dcs07/afeinber/data/personal/jzhan/CPEL2/ROIcut.R "$ABSOLUTE_MAIN_DIR"


# Loop through chromosomes 1-19
for chr in {1..19}; do
    # Create chromosome subdirectory
    CHR_DIR="chr_${chr}"
    mkdir -p "$CHR_DIR"
    
    echo "Processing chromosome $chr"
    echo "Created subdirectory: $CHR_DIR"
    
    # Run Julia script with chromosome number and output directory

    julia --project=/dcs04/feinberg/data/shared/CPELTDM -p 8 /dcs07/afeinber/data/personal/jzhan/CPEL2/main.jl "$CHR_DIR" "$chr" "$ABSOLUTE_MAIN_DIR"
    
    
     
    echo "Completed processing chromosome $chr"
done

echo "All chromosomes processed successfully!"

