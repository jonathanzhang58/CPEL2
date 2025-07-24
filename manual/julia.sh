#!/bin/bash
#SBATCH --mail-user=jzhan367@jh.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cpel_parallelcut 
#SBATCH --chdir=.
#SBATCH --partition=cee,cegs2,shared
#SBATCH --mem=50G
#SBATCH --time=400:00:00
#SBATCH --output=./manual/logfiles/julia_%j.out
#SBATCH --error=./manual/logfiles/julia_%j.err
#SBATCH --cpus-per-task=8
# Read command line arguments

ABSOLUTE_MAIN_DIR=${1} 

cd "$ABSOLUTE_MAIN_DIR" || exit

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
      /dcs07/afeinber/data/personal/jzhan/CPEL2/scripts/sh_scripts/chr_anal.sh 8 "$CHR_DIR" "$chr" "$ABSOLUTE_MAIN_DIR"

   echo "Began processing chromosome $chr"
done

echo "All chromosomes started processed successfully!"
