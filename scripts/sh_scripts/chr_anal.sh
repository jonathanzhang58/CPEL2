#!/bin/bash
# No email notifications
#SBATCH --mail-type=FAIL,END
#SBATCH --chdir=.
#SBATCH --partition=cee,cegs2,shared
#SBATCH --mem=40G
#SBATCH --time=400:00:00
#SBATCH --cpus-per-task=8


# Read command line arguments
NUM_CORES=${1:-8}
OUTDIR=${2}
CHR_NUM=${3}
INDIR=${4}
FA=${5}
BAMS1=${6}
BAMS2=${7}

julia --project=/dcs04/feinberg/data/shared/CPELTDM -p "$NUM_CORES" /dcs07/afeinber/data/personal/jzhan/CPEL2/scripts/main.jl "$OUTDIR" "$CHR_NUM" "$INDIR" "$FA" "$BAMS1" "$BAMS2"  