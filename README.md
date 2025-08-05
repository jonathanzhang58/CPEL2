# CpelTdm Pipeline

A comprehensive bioinformatics pipeline for identifying differentially methylated regions (DMRs) from bisulfite sequencing data using the CpelTDM package

## Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Pipeline Architecture](#pipeline-architecture)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Technical Details](#technical-details)
- [Contributing](#contributing)

## Overview

The CPEL2 pipeline performs analysis of DNA methylation patterns on a subset of the genome, in order to increase power. The package can detect DMRs based on mean methylation level (MML), normalized
methylation entropy (NME), and the coefficient of methylation divergence between probability distributions
(PDMs) of methylation patterns. 

For further information, see the original CpelTdm paper: https://www.biorxiv.org/content/10.1101/2020.10.17.343020v1.full.pdf

Or repository: 
https://github.com/jordiabante/CpelTdm.jl

The pipeline integrates:

- **ROI Creation**: Generation of regions of interest based on genomic landmarks (promoters, CpG shores, CpG islands)
- **CPEL Analysis**: Parallelized chromosome-specific differential methylation analysis for increased spped. 
- **Data Consolidation**: Genome-wide integration of results
- **Statistical Filtering**: Adjusted P-value based filtering of significant DMRs
- **Further Analysis**: All signficant intermediate files are saved for further analysis purposes. 

### Expected Usage
The pipeline is meant for usage on the JHPCE cluster, which uses the SLURM scheduler. It can be easily adapted to other usages. 
## Quick Start

### 1. Clone and Setup
```bash
git clone <repository-url>
cd CPEL2
```

### 2. Configure Data Paths
Edit `datafile.txt` with your data locations:
```
/path/to/your/bsseq/data.RDS
/path/to/group1/bam/files
/path/to/group2/bam/files
/path/to/reference/genome.fa
```
### 3. Set Working Directory
Edit `datafile.txt` with your data locations:
```
cd ./CPELPipeline
```

### 4. Run Pipeline
```bash
# Basic analysis (promoters, 250bp tiles, p-value 0.1)
sbatch main.sh

# Custom parameters
sbatch main.sh shores 500 0.05
```

## Installation

### Prerequisites

#### System Requirements
- **SLURM**: Job scheduling system
- **Linux/Unix**: Operating system
- **Memory**: 50GB RAM recommended
- **Storage**: Sufficient space for intermediate and final results
- **CPU**: Multi-core system for parallel processing

#### Software Dependencies
- **R** (≥4.0.0) with Bioconductor packages:
  - `bsseq` - Bisulfite sequencing analysis
  - `GenomicRanges` - Genomic data manipulation
  - `AnnotationDbi` - Annotation database interface
  - `rtracklayer` - Genomic data import/export
  - `Mus.musculus` - Mouse genome annotations
  - `TxDb.Mmusculus.UCSC.mm10.knownGene` - Mouse gene annotations
- **Julia** (≥1.6.0) - For chromosome-specific analysis
- **Bash** - Shell scripting

### Installation Steps

1. **Clone Repository**
   ```bash
   git clone <repository-url>
   cd CPEL2
   ```

2. **Install R Dependencies**
   ```r
   # Install BiocManager if not already installed
   if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
   
   # Install required Bioconductor packages
   BiocManager::install(c("bsseq", "GenomicRanges", "AnnotationDbi", 
                         "rtracklayer", "Mus.musculus", 
                         "TxDb.Mmusculus.UCSC.mm10.knownGene"))
   ```

3. **Verify Installation**
   ```bash
   # Test R dependencies
   Rscript -e "library(bsseq); library(GenomicRanges); cat('R dependencies OK\n')"
   
   # Test Julia
   julia --version
   ```

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
- Files should be named consistently (e.g., `sample1.bam`, `sample2.bam`)

**Reference Genome**:
- FASTA format
- Should match the genome used for alignment
- Index files (`.fai`) should be present

### Reference Files

Place required reference files in the `data/` directory:

```
CPEL2/
├── data/
│   ├── CpG_islands_mm10.txt    # CpG islands annotation (required)
│   └── other_reference_files/  # Additional reference files
```

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

### Job Monitoring

```bash
# Check job status
squeue -u $USER

# Monitor specific job
squeue -j <job_id>

# View job logs
tail -f output/logfiles/cpel_main_pipeline*.out
tail -f output/logfiles/cpel_analyze_consolidated*.out
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

### Pipeline Flow

```
1. ROI Creation
   ├── Read datafile.txt configuration
   ├── Extract genomic landmarks (promoters/shores/islands)
   ├── Create tiled regions of interest
   └── Export BED files

2. Chromosome Analysis (Parallel)
   ├── Submit jobs for chromosomes 1-19
   ├── Process each chromosome independently
   └── Generate intermediate results

3. Data Consolidation
   ├── Wait for all chromosome jobs to complete
   ├── Combine results from all chromosomes
   └── Apply statistical filtering

4. Final Output
   ├── Filter by p-value threshold
   ├── Export filtered DMRs to CSV
   └── Generate summary statistics
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
    ├── chr_1/, chr_2/, ..., chr_19/    # Chromosome-specific results
    ├── RDS_files<timestamp>/           # Individual RDS files
    ├── concatenated_output<timestamp>/ # Combined results
    │   └── filtered_output<timestamp>/ # Final filtered results
    └── logfiles/                       # Job execution logs
```

### File Types

#### BED Files
- **ROI BED files**: Regions of interest in BED format
- **Filtered BED files**: Final DMR regions

#### RDS Files
- **Individual RDS**: Chromosome-specific analysis results
- **Concatenated RDS**: Genome-wide combined results

#### CSV Files
- **Filtered CSV**: Final DMR results in human-readable format
- **Summary CSV**: Statistical summary of results

#### Log Files
- **Main pipeline logs**: Overall pipeline execution
- **Analysis logs**: Individual analysis job logs
- **Error logs**: Error messages and debugging information

### Output Interpretation

#### DMR Results (CSV)
Columns include:
- `seqnames`: Chromosome name
- `start`, `end`: Genomic coordinates
- `width`: Region width in base pairs
- `score`: Methylation score
- `p.val`: Raw p-value
- `p.adj`: Adjusted p-value (FDR-corrected)
- Additional metadata columns

#### Summary Statistics
- Total regions analyzed
- Number of significant DMRs
- Coverage statistics
- Quality metrics

## Troubleshooting

### Common Issues

#### 1. Missing Data Files
```
Error: datafile.txt not found
```
**Solution**: Ensure `datafile.txt` exists in the project root and contains valid paths.

#### 2. Missing Reference Files
```
Error: CpG islands file not found
```
**Solution**: Place `CpG_islands_mm10.txt` in the `data/` directory.

#### 3. R Package Errors
```
Error: package 'bsseq' not available
```
**Solution**: Install required Bioconductor packages:
```r
BiocManager::install(c("bsseq", "GenomicRanges", "AnnotationDbi"))
```

#### 4. Memory Issues
```
Error: cannot allocate vector of size X
```
**Solution**: Increase memory allocation in SLURM script or reduce batch size.

#### 5. Job Failures
```
Error: SLURM job failed
```
**Solution**: Check job logs in `output/logfiles/` for specific error messages.

### Debugging

#### Enable Verbose Logging
```bash
# Add debug output to scripts
export DEBUG=1
sbatch main.sh
```

#### Check Job Status
```bash
# View recent jobs
sacct -u $USER --starttime=2024-01-01

# Check job details
scontrol show job <job_id>
```

#### Validate Input Data
```r
# Test BSseq object
library(bsseq)
bs_obj <- readRDS("path/to/bsseq.RDS")
print(bs_obj)
```

### Performance Optimization

#### Memory Usage
- **Large datasets**: Increase `--mem` in SLURM scripts
- **Many samples**: Consider subsetting data for testing
- **High coverage**: Adjust coverage thresholds in scripts

#### Runtime Optimization
- **Parallel processing**: Ensure sufficient CPU cores
- **I/O optimization**: Use fast storage for temporary files
- **Batch size**: Adjust based on available resources

## Technical Details

### Algorithm Overview

1. **ROI Generation**
   - Extract genomic landmarks (promoters: ±2kb from TSS, shores: ±2kb from CpG islands)
   - Create fixed-width tiles within landmark regions
   - Filter tiles by CpG density and coverage requirements

2. **CPEL Analysis**
   - Calculate methylation levels for each ROI
   - Perform statistical testing between groups
   - Apply multiple testing correction

3. **Quality Control**
   - Coverage filtering (≥10x in ≥5 samples per group)
   - CpG density filtering (≥3 CpGs per tile)
   - Distance filtering (≤1kb from genes for promoters)

### Statistical Methods

- **Differential methylation**: Welch's t-test or similar
- **Multiple testing correction**: Benjamini-Hochberg FDR
- **Coverage filtering**: Minimum coverage thresholds per group
- **Quality metrics**: CpG density, distance to features

### Data Formats

#### Input Formats
- **BSseq RDS**: R serialized BSseq object
- **BAM**: Binary alignment files (sorted, indexed)
- **FASTA**: Reference genome sequence
- **BED**: Genomic annotations (CpG islands)

#### Output Formats
- **BED**: Genomic coordinates and scores
- **CSV**: Tabular results with statistics
- **RDS**: R objects for further analysis

### Version Information

- **Pipeline Version**: 2.0
- **R Version**: ≥4.0.0
- **Bioconductor**: ≥3.15
- **Julia Version**: ≥1.6.0
- **Last Updated**: 2024

## Contributing

### Development Setup

1. **Fork Repository**
   ```bash
   git clone <your-fork-url>
   cd CPEL2
   ```

2. **Create Feature Branch**
   ```bash
   git checkout -b feature/new-feature
   ```

3. **Make Changes**
   - Follow existing code style
   - Add tests for new functionality
   - Update documentation

4. **Submit Pull Request**
   - Include detailed description
   - Reference related issues
   - Ensure all tests pass

### Code Style

- **Bash**: Use consistent indentation (2 spaces)
- **R**: Follow tidyverse style guide
- **Julia**: Follow Julia style conventions
- **Documentation**: Update README for new features

### Testing

```bash
# Run basic tests
./test_pipeline.sh

# Validate configuration
Rscript scripts/validate_config.R

# Check output format
Rscript scripts/validate_output.R
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this pipeline in your research, please cite:

```
CPEL2 Pipeline: A Comprehensive Tool for Differential Methylation Analysis
Authors: [Your Name]
Year: 2024
```

## Support

For questions and support:

- **Issues**: Create an issue on GitHub
- **Documentation**: Check this README and inline comments
- **Email**: [Your Email]

---

**Note**: This pipeline is designed for research use. Always validate results and consider biological context when interpreting DMRs.
