# CpelTdm Pipeline

A comprehensive bioinformatics pipeline for identifying differentially methylated regions (DMRs) from bisulfite sequencing data using the CpelTDM package

## Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
- [Dependencies](#dependencies)
- [Configuration](#configuration)
- [Usage](#usage)
- [Architecture](#pipeline-architecture)
- [Output Files](#output-files)


## Overview

The CPEL2 pipeline performs analysis of DNA methylation patterns on a subset of the genome, in order to increase power. The package can detect DMRs based on mean methylation level (MML), normalized
methylation entropy (NME), and the coefficient of methylation divergence between probability distributions
(PDMs) of methylation patterns. 

The pipeline accomplishes this by first filtering the genome such that only those regions around specific genomic landmarks remain (promoters, CpG shores, CpG islands). The script then breaks up these regions into tiles (default 250bp) and further filters by removing tiles with insufficient coverage across samples or insufficient CpG site density. Then CPEL is used on the remaining tiles (exported as .bed) to compute various scores that are used to find DMRs. The script then filters all identified DMRs to only those with sufficiently high adjusted p-value (0.1 default).  

For further information, see the original CpelTdm paper: https://www.biorxiv.org/content/10.1101/2020.10.17.343020v1.full.pdf

Or repository: 
https://github.com/jordiabante/CpelTdm.jl


### Expected Usage
The pipeline is meant for usage on the JHPCE cluster, which uses the SLURM scheduler. It can be easily adapted to other usages. 
## Quick Start

### 1. Install Prerequisites
Install Julia 1.3.1 by going to: 
```https://julialang.org/downloads/oldreleases/```

and downloading the appropriate `tar.gz` file and moving that to your home directory in the JHPCE using your method choice.



Replace text of file at 

```
/yourhome/directory/.julia/environments/v1.3/Mainfest.toml
```

With the Manifest.toml in the repository


Then run (in Julia): 

```
Pkg.resolve()
Pkg.instantiate()
```


### 2. Clone and Setup
```bash
git clone https://github.com/jonathanzhang58/CpelTdm_pipeline
cd CpelTdm_pipeline  
```

### 3. Configure Data Paths
Edit `datafile.txt` with your data locations:
```
/path/to/your/bsseq/data.RDS
/path/to/group1/bam/files
/path/to/group2/bam/files
/path/to/reference/genome.fa
```
### 4. Set Working Directory
Edit `datafile.txt` with your data locations:
```
cd ./CPELPipeline
```

### 5. Run Pipeline
```bash
# Basic analysis (promoters, 250bp tiles, p-value 0.1)
sbatch main.sh

# Custom parameters
sbatch main.sh shores 500 0.05
```

## Dependencies

#### Software Dependencies
- **R** (≥3.9.0) with Bioconductor packages:
  - `bsseq` - Bisulfite sequencing analysis
  - `GenomicRanges` - Genomic data manipulation
  - `AnnotationDbi` - Annotation database interface
  - `rtracklayer` - Genomic data import/export
  - `Mus.musculus` - Mouse genome annotations
  - `TxDb.Mmusculus.UCSC.mm10.knownGene` - Mouse gene annotations
- **Julia** (1.3.1) - For chromosome-specific analysis. MUST BE THIS SPECIFIC VERSION
- **Bash** - Shell scripting

## Configuration

### Data File Configuration

The `datafile.txt` file contains paths to your input data:

```
Line 1: /path/to/bsseq/data.RDS          # BSseq object path
Line 2: /path/to/group1/bam/files        # Group 1 BAM files directory
Line 3: /path/to/group2/bam/files        # Group 2 BAM files directory  
Line 4: /path/to/reference/genome.fa     # Reference genome FASTA
```

#### Data Format Requirements

**BSseq Object**:
- Must be an RDS file containing a BSseq object
- Should include coverage and methylation data for all samples
- Sample names should match those in BAM directories

**BAM Files**:
- Sorted and indexed BAM files
- One directory per experimental group
- Index files (`.bai`) should be present


**Reference Genome**:
- FASTA format
- Should match the genome used for alignment
- Index files (`.fai`) should be present


## Usage

### Basic Usage

```bash
sbatch main.sh [LANDMARK] [WIDTH] [PVAL]
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `LANDMARK` | string | `promoters` | Genomic region type: `promoters`, `shores`, `islands` |
| `WIDTH` | integer | `250` | Tile width in base pairs |
| `PVAL` | float | `0.1` | P-value threshold for filtering |

### Usage Examples

```bash
# Default analysis (promoters, 250bp tiles, p-value 0.1)
sbatch main.sh

# CpG shores analysis with custom parameters
sbatch main.sh shores 500 0.05

# CpG islands analysis with strict filtering
sbatch main.sh islands 250 0.01

# Promoter analysis with large tiles
sbatch main.sh promoters 1000 0.1
```


## Pipeline Architecture

### Directory Structure

```
CPEL2/
├── main.sh                          # Main pipeline entry point
├── datafile.txt                     # Data configuration file
├── README.md                        # This documentation
├── CpG_islands_mm10.txt            # CpG islands reference file
├── scripts/
│   ├── consolidated_analysis.R      # Main analysis script
│   ├── ROIcut.R                     # ROI processing script
│   ├── main.jl                      # Julia analysis script
│   ├── partitioning/
│   │   ├── partitioning.R           # ROI creation script
│   │   └── filtering.R              # Landmark extraction functions
│   └── sh_scripts/
│       ├── analyze_consolidated.sh  # Analysis job submission
│       └── chr_anal.sh             # Chromosome-specific analysis
├── output/
│   └── logfiles/                   # Job execution logs
└── deprecated_files/               # Legacy scripts (not used)
```

### Script Descriptions

#### Main Scripts

- **`main.sh`**: Orchestrates the complete pipeline workflow
- **`consolidated_analysis.R`**: Performs analysis, merging, and filtering in a single R session
- **`partitioning.R`**: Creates regions of interest based on genomic landmarks
- **`ROIcut.R`**: Processes and prepares ROIs for analysis

#### Supporting Scripts

- **`filtering.R`**: Functions for extracting different genomic landmarks
- **`chr_anal.sh`**: Submits chromosome-specific analysis jobs
- **`analyze_consolidated.sh`**: Submits the final consolidation job

## Output Files

### Directory Structure

Results are organized in timestamped directories:

```
output/
└── unfiltered_cpel_<landmark>_<width>_<timestamp>/
    ├── *.bed                           # Initial ROI BED files
    ├── chr_1/, chr_2/, ..., chr_19/    # Chromosome-specific results and logfiles
    ├── RDS_files<timestamp>/           # Individual RDS files
    ├── concatenated_output<timestamp>/ # Combined results
    │   └── filtered_output<timestamp>/ # Final filtered results
    └── logfiles/                       # Job execution logs
```

### Output Interpretation

#### DMR Results (CSV)
Columns include:
- `seqnames`: Chromosome name
- `start`, `end`: Genomic coordinates
- `width`: Region width in base pairs
- `score`: Methylation score
- `p.val`: Raw p-value
- `p.adj`: Adjusted p-value
- Additional metadata columns

#### Summary Statistics
- Total regions analyzed
- Number of significant DMRs
- Coverage statistics
- Quality metrics