#!/bin/bash
#SBATCH --mail-user=jzhan367@jh.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=cpel_analyze
#SBATCH --chdir=.
#SBATCH --partition=shared
#SBATCH --mem=5G
#SBATCH --time=400:00:00
#SBATCH --output=./output/logfiles/cpel_analyze%j.out
#SBATCH --error=./output/logfiles/cpel_analyze%j.err
#SBATCH --cpus-per-task=8
module load conda_R
set -euo pipefail

# Read command line arguments
in_dir=${1}
p_val=${2}

echo "Analyzing $in_dir with p-value $p_val" 
# Get current date in yy_mm_dd_hh_mm_ss format
DATE=$(date +%y_%m_%d_%H_%M_%S)

# Create main output directory
outdir="${in_dir}/RDS_files${DATE}"
mkdir -p "$outdir"


# Loop through each subdirectory
for subdir in "$in_dir"/*/; do
    # Remove trailing slash and get the directory name
    dir_name=$(basename "${subdir%/}")
    sub_outdir="${outdir}/${dir_name}_out"
    mkdir -p "$sub_outdir"

    # For each bedGraph file in the subdirectory
    for bedgraph_file in "${subdir%/}"/*.bedGraph; do
        if [[ -f "$bedgraph_file" ]]; then
            echo "Processing $bedgraph_file -> $sub_outdir"
            Rscript /dcs07/afeinber/data/personal/jzhan/CPEL2/scripts/Analyze.R "$bedgraph_file" "$sub_outdir"
        fi
    done
done

concat_outdir="${in_dir}/concatenated_output${DATE}"
mkdir -p "$concat_outdir"

Rscript /dcs07/afeinber/data/personal/jzhan/CPEL2/scripts/merge.R "$outdir" "$concat_outdir"

sbatch /dcs07/afeinber/data/personal/jzhan/CPEL2/scripts/sh_scripts/filter.sh "$concat_outdir" "$p_val"