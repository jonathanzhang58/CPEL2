#!/bin/bash
#SBATCH --mail-user=jzhan367@jh.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cpel_analyze
#SBATCH --chdir=.
#SBATCH --partition=shared
#SBATCH --mem=200G
#SBATCH --time=400:00:00
#SBATCH --output=cpel_analyze.out
#SBATCH --error=cpel_analyze.err
#SBATCH --cpus-per-task=8
module load conda_R
set -euo pipefail



# Get current date in yy_mm_dd_hh_mm_ss format
DATE=$(date +%y_%m_%d_%H_%M_%S)

# Create main output directory
outdir="analyzed_output${DATE}"
mkdir -p "$outdir"


# Loop through each subdirectory in cpel_out/target/
for subdir in cpel_out/target/*/; do
    # Remove trailing slash and get the directory name
    dir_name=$(basename "${subdir%/}")
    sub_outdir="${outdir}/${dir_name}_out"
    mkdir -p "$sub_outdir"

    # For each bedGraph file in the subdirectory
    for bedgraph_file in "${subdir%/}"/*.bedGraph; do
        if [[ -f "$bedgraph_file" ]]; then
            echo "Processing $bedgraph_file -> $sub_outdir"
            Rscript Analyze.R "$bedgraph_file" "$sub_outdir"
        fi
    done
done

mkdir -p "concatenated_output"
Rscript statistics.R "$outdir" "concatenated_output"