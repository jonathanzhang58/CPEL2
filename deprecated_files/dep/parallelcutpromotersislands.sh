#!/bin/bash
#SBATCH --mail-user=jzhan367@jh.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cpel_parallelcut 
#SBATCH --chdir=.
#SBATCH --partition=cee,cegs2,shared
#SBATCH --mem=50G
#SBATCH --time=400:00:00
#SBATCH --output=./output/logfiles/islandsparallelcut%j.out
#SBATCH --error=./output/logfiles/islandsparallelcut_%j.err
#SBATCH --cpus-per-task=8
# Read command line arguments

WIDTH=${1:-250}  # Default to 250 if not provided

module load conda_R
cd /dcs07/afeinber/data/personal/jzhan/CPEL2/output || exit
# Get current date in yy_mm_dd_hh_mm_ss format
DATE=$(date +%y_%m_%d_%H_%M_%S)

# Create main output directory
MAIN_DIR="islands_cpel_old_young_${WIDTH}_${DATE}"
mkdir -p "$MAIN_DIR"

echo "Created main directory: $MAIN_DIR"

# Store the absolute path before changing directories
ABSOLUTE_MAIN_DIR="$(pwd)/$MAIN_DIR"

cd "$MAIN_DIR" || exit

echo "Creating ROIs"
Rscript /dcs07/afeinber/data/personal/jzhan/CPEL2/250bpislands.R "$WIDTH"
# Check if any .bed files exist in the current directory
if ! ls *.bed 1> /dev/null 2>&1; then
    echo "Error: No .bed files found in $MAIN_DIR"
    exit 1
fi

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
      -J "cut_shores_cpel_chr_anal_${chr}" \
      -o "$CHR_DIR/logfiles/test_cpel_%j.out" \
      -e "$CHR_DIR/logfiles/test_cpel_%j.err" \
      /dcs07/afeinber/data/personal/jzhan/CPEL2/chr_anal.sh 8 "$CHR_DIR" "$chr" "$ABSOLUTE_MAIN_DIR"

   echo "Began processing chromosome $chr"
done

echo "All chromosomes started processed successfully!"
