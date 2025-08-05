#!/bin/bash
#SBATCH --mail-user=your_email@jh.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cpel_main_pipeline
#SBATCH --chdir=.
#SBATCH --partition=shared
#SBATCH --mem=50G
#SBATCH --time=400:00:00
#SBATCH --output=output/logfiles/cpel_main_pipeline%j.out
#SBATCH --error=output/logfiles/cpel_main_pipeline%j.err
#SBATCH --cpus-per-task=8

# CPEL Main Pipeline Script
# =========================
# This script orchestrates the complete CPEL pipeline:
# 1. Creates Regions of Interest (ROIs) based on user-specified genomic landmarks
# 2. Runs CPEL analysis on each chromosome to find differentially methylated regions (DMRs)
# 3. Combines outputs from all chromosomes into a single genome-wide file
# 4. Filters results based on p-value threshold
#
# Prerequisites:
# - datafile.txt must exist with the following format:
#   Line 1: path/to/bsseq
#   Line 2: path/to/group1/bams
#   Line 3: path/to/group2/bams
#   Line 4: path/to/fasta
#
# Usage: sbatch main.sh [LANDMARK] [WIDTH] [PVAL]
#   LANDMARK: Type of genomic region (promoters/shores/islands, default: promoters)
#   WIDTH:    Tile width in base pairs (default: 250)
#   PVAL:     P-value threshold for filtering (default: 0.1)
#
# Example: sbatch main.sh shores 500 0.05

LANDMARK=${1:-PROMOTERS}
WIDTH=${2:-250}  # Default to 250 if not provided
PVAL=${3:-0.1}   # Default p-value threshold to 0.1 if not provided




# Get the directory where this script is located (project root)

# SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_ROOT="$(pwd)"

echo "Project root: $PROJECT_ROOT"

# Read configuration from datafile.txt (relative to project root)
DATAFILE="$PROJECT_ROOT/datafile.txt"
if [[ ! -f "$DATAFILE" ]]; then
    echo "Error: datafile.txt not found at $DATAFILE"
    echo "Please create datafile.txt in the project root directory"
    exit 1
fi

echo "Reading configuration from: $DATAFILE"
# Read each line from datafile.txt into variables
BSSEQ_PATH=$(sed -n '1p' "$DATAFILE")
GROUP1_BAMS_PATH=$(sed -n '2p' "$DATAFILE")
GROUP2_BAMS_PATH=$(sed -n '3p' "$DATAFILE")
FASTA_PATH=$(sed -n '4p' "$DATAFILE")

# Validate that paths are not empty
if [[ -z "$BSSEQ_PATH" || -z "$GROUP1_BAMS_PATH" || -z "$GROUP2_BAMS_PATH" || -z "$FASTA_PATH" ]]; then
    echo "Error: One or more paths in datafile.txt are empty"
    echo "Expected format:"
    echo "Line 1: path/to/bsseq"
    echo "Line 2: path/to/group1/bams"
    echo "Line 3: path/to/group2/bams"
    echo "Line 4: path/to/fasta"
    exit 1
fi

# Convert LANDMARK to lowercase
LANDMARK=$(echo "$LANDMARK" | tr '[:upper:]' '[:lower:]')

module load conda_R

# Ensure output and logfiles directories exist
mkdir -p "$PROJECT_ROOT/output/logfiles"

cd "$PROJECT_ROOT/output" || exit
# Get current date in yy_mm_dd_hh_mm_ss format
DATE=$(date +%y_%m_%d_%H_%M_%S)

# Create main output directory
MAIN_DIR="unfiltered_cpel_old_young_${WIDTH}_${DATE}"
mkdir -p "$MAIN_DIR"

echo "Created main directory: $MAIN_DIR"

# Store the absolute path before changing directories
ABSOLUTE_MAIN_DIR="$(pwd)/$MAIN_DIR"

cd "$MAIN_DIR" || exit

echo "Creating ROIs"

Rscript "$PROJECT_ROOT/scripts/partitioning/partitioning.R" "$WIDTH" "$LANDMARK" "$BSSEQ_PATH"
# Check if any .bed files exist in the current directory
if ! ls *.bed 1> /dev/null 2>&1; then
    echo "Error: No .bed files found in $MAIN_DIR"
    exit 1
fi


# Run ROIcut.R with the absolute main directory as argument
echo "Running ROIcut.R with directory: $ABSOLUTE_MAIN_DIR"
Rscript "$PROJECT_ROOT/scripts/ROIcut.R" "$ABSOLUTE_MAIN_DIR"
job_ids=""
# Loop through chromosomes 1-19
for chr in {1..19}; do
   # Create chromosome subdirectory
   CHR_DIR="chr_${chr}"
   mkdir -p "$CHR_DIR"
   echo "Created subdirectory: $CHR_DIR"
   mkdir -p "$CHR_DIR/logfiles"
   echo "Created logfiles subdirectory: $CHR_DIR/logfiles"

   # Run Julia script with chromosome number and output directory
  jid=$(sbatch --parsable \
      -J "chr_anal_${chr}" \
      -o "$CHR_DIR/logfiles/test_cpel_%j.out" \
      -e "$CHR_DIR/logfiles/test_cpel_%j.err" \
      "$PROJECT_ROOT/scripts/sh_scripts/chr_anal.sh" 8 "$CHR_DIR" "$chr" "$ABSOLUTE_MAIN_DIR" "$FASTA_PATH" "$GROUP1_BAMS_PATH" "$GROUP2_BAMS_PATH" "$PROJECT_ROOT")
  job_ids="${job_ids}:${jid}"

done

echo "All chromosomes started processed successfully!"

# Remove leading colon
job_ids="${job_ids#:}"

# Submit consolidated analysis after all chromosome jobs
echo "Submitting consolidated analysis with dependency on jobs: ${job_ids}"
analysis_job_id=$(sbatch --parsable --dependency=afterok:"${job_ids}" \
    "$PROJECT_ROOT/scripts/sh_scripts/analyze_consolidated.sh" "$ABSOLUTE_MAIN_DIR" "$PVAL" "$PROJECT_ROOT")

echo "All chromosome jobs submitted successfully!"
echo "Analysis job ID: $analysis_job_id"
echo ""
echo "Pipeline Summary:"
echo "=================="
echo "Landmark type: $LANDMARK"
echo "Tile width: ${WIDTH}bp"
echo "P-value threshold: $PVAL"
echo "Output directory: $ABSOLUTE_MAIN_DIR"
echo ""
echo "Data Paths:"
echo "----------"
echo "BSseq data: $BSSEQ_PATH"
echo "Group 1 BAMs: $GROUP1_BAMS_PATH"
echo "Group 2 BAMs: $GROUP2_BAMS_PATH"
echo "FASTA reference: $FASTA_PATH"
echo ""
echo "The pipeline will:"
echo "1. Create ROIs based on $LANDMARK with ${WIDTH}bp tiles"
echo "2. Run CPEL analysis on each chromosome (jobs: ${job_ids})"
echo "3. Combine outputs and filter by p-value $PVAL (job: $analysis_job_id)"
echo ""
echo "Monitor progress with: squeue -u \$USER"
echo "Final results will be in: $ABSOLUTE_MAIN_DIR/concatenated_output*/filtered_output*/"