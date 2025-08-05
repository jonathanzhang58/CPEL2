#!/bin/bash
#SBATCH --mail-user=jzhan367@jh.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cpel_stats
#SBATCH --chdir=.
#SBATCH --partition=shared
#SBATCH --mem=200G
#SBATCH --time=400:00:00
#SBATCH --output=cpel_stats.out
#SBATCH --error=cpel_stats.err
#SBATCH --cpus-per-task=8
module load conda_R
set -euo pipefail

# Read command line arguments
in_dir=${1}

# Get current date in yy_mm_dd_hh_mm_ss format
DATE=$(date +%y_%m_%d_%H_%M_%S)

# Create main output directory
outdir="${in_dir}/analyzed_concatenated_output${DATE}"
mkdir -p "$outdir"


Rscript statistics.R "$in_dir" "$outdir"

