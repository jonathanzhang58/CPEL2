#!/bin/bash
#SBATCH --mail-user=jzhan367@jh.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cpel_analyze_consolidated
#SBATCH --chdir=.
#SBATCH --partition=shared
#SBATCH --mem=50G
#SBATCH --time=400:00:00
#SBATCH --output=output/logfiles/cpel_analyze_consolidated%j.out
#SBATCH --error=output/logfiles/cpel_analyze_consolidated%j.err
#SBATCH --cpus-per-task=8
module load conda_R
set -euo pipefail

# Read command line arguments
in_dir=${1}
p_val=${2}
project_root=${3:-$(pwd)}

echo "Running consolidated analysis on $in_dir with p-value $p_val" 

# Run the consolidated script that performs all operations in memory
# This replaces lines 27-49 from the original analyze.sh
Rscript "$project_root/scripts/consolidated_analysis.R" "$in_dir" "$p_val"

echo "Consolidated analysis complete!"