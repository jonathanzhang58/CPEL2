#!/bin/bash
#SBATCH --mail-user=jzhan367@jh.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cpel_parallel
#SBATCH --chdir=.
#SBATCH --partition=cee,cegs2,shared
#SBATCH --mem=200G
#SBATCH --time=400:00:00
#SBATCH --output=./output/logfiles/parallel%j.out
#SBATCH --error=./output/logfiles/parallel_%j.err
#SBATCH --cpus-per-task=8

module load conda_R
cd /dcs07/afeinber/data/personal/jzhan/CPEL2/output || exit
# Get current date in yy_mm_dd_hh_mm_ss format
DATE=$(date +%y_%m_%d_%H_%M_%S)

# Create main output directory
MAIN_DIR="cpel_old_young_${DATE}"
mkdir -p "$MAIN_DIR"

echo "Created main directory: $MAIN_DIR"

# Store the absolute path before changing directories
ABSOLUTE_MAIN_DIR="$(pwd)/$MAIN_DIR"

cd "$MAIN_DIR" || exit

echo "Creating ROIs"
Rscript /dcs07/afeinber/data/personal/jzhan/CPEL2/shores.R

# Run ROIcut.R with the absolute main directory as argument
echo "Running ROIcut.R with directory: $ABSOLUTE_MAIN_DIR"
Rscript /dcs07/afeinber/data/personal/jzhan/CPEL2/ROIcut.R "$ABSOLUTE_MAIN_DIR"

# Loop through chromosomes 1-19
for chr in {1..19}; do
   # Create chromosome subdirectory
   CHR_DIR="chr_${chr}"
   mkdir -p "$CHR_DIR"
   echo "Created subdirectory: $CHR_DIR"
   mkdir -p "$CHR_DIR/logfiles"
   echo "Created logfiles subdirectory: $CHR_DIR/logfiles"

   # Run Julia script with chromosome number and output directory

   sbatch \
      -J "cpel_chr_anal_${chr}" \
      -o "$CHR_DIR/logfiles/test_cpel_%j.out" \
      -e "$CHR_DIR/logfiles/test_cpel_%j.err" \
      /dcs07/afeinber/data/personal/jzhan/CPEL2/chr_anal.sh 8 "$CHR_DIR" "$chr" "$ABSOLUTE_MAIN_DIR"

   echo "Began processing chromosome $chr"
done

echo "All chromosomes started processed successfully!"
