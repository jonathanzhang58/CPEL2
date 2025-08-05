#!/usr/bin/env bash
#SBATCH --mail-user=jzhan367@jh.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=test_cpel
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=8
#SBATCH --partition=cee,cegs2,shared
#SBATCH --time=96:00:00
#SBATCH --output=./logfiles/test_cpel_%j.out
#SBATCH --error=./logfiles/test_cpel_%j.err

base_dir="/dcs07/afeinber/data/personal/jzhan/CPEL2/jonathan_test"

julia --project=/dcs04/feinberg/data/shared/CPELTDM -p 8 ${base_dir}/jonathan_test.jl ${base_dir}/test_output chr19 ${base_dir}
