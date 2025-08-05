#!/usr/bin/env bash
#SBATCH --mail-user=jzhan367@jh.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=sort_reheader_bams
#SBATCH --mem=80G
#SBATCH --time=96:00:00
#SBATCH --output=sort_reheader_bams.out
#SBATCH --error=sort_reheader_bams.err


module load samtools

fl=$1


samtools view -H $fl | sed -E 's/SN:\bchr([0-9]+)\b/SN:\1/g' | samtools reheader - $fl > "${fl}_noCHR.bam"
samtools index "${fl}_noCHR.bam"

