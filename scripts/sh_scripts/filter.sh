#!/bin/bash
#SBATCH --mail-user=jzhan367@jh.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cpel_filter
#SBATCH --chdir=.
#SBATCH --partition=shared
#SBATCH --mem=50G
#SBATCH --time=400:00:00
#SBATCH --output=./output/logfiles/cpel_filter%j.out
#SBATCH --error=./output/logfiles/cpel_filter%j.err
#SBATCH --cpus-per-task=8
module load conda_R
set -euo pipefail

# Read command line arguments
in_dir=${1}
p_val=${2:-0.1}

echo "Filtering $in_dir with p-value $p_val" 

# Get current date in yy_mm_dd_hh_mm_ss format
DATE=$(date +%y_%m_%d_%H_%M_%S)

# Create main output directory
outdir="${in_dir}/filtered_output${DATE}"
mkdir -p "$outdir"


Rscript /dcs07/afeinber/data/personal/jzhan/CPEL2/scripts/filter_sum.R "$in_dir" "$outdir" "$p_val"