using Distributed
@everywhere using CpelTdm

# Input arguments will be the desired output directory and the chromosome number (1-22). 
outdir = ARGS[1]
chr = ARGS[2]
indir = ARGS[3]
fa = ARGS[4]
bams1 = ARGS[5]
bams2 = ARGS[6]

println("Output directory: ", outdir)
println("Chromosome: ", chr)
println("Input directory: ", indir)

# Genome file
fa = ARGS[4]

# Give it the file containing your ROIs.
bed = "$(indir)/chr/myROIS_chr$(chr).bed"
println("BED file: ", bed)

# Check if BED file exists
if !isfile(bed)
    println("ERROR: BED file does not exist: ", bed)
    exit(1)
end

# Put the paths to sorted.bam files for each group in two lists.
# Get list of full paths:
bam_dir = "/dcs05/feinberg/data/personal/ocamacho/WGBS_aging_epidermis/informME/coord_sorted_bams"

# Check if BAM directory exists
if !isdir(bam_dir)
    println("ERROR: BAM directory does not exist: ", bam_dir)
    exit(1)
end

# Old
bams1 = readdir(bam_dir)
# Young 
bams2 = readdir(bam_dir)



if isempty(bams1) || isempty(bams2)
    println("ERROR: No BAM files found in directory")
    println("Old BAM files found: ", length(bams1))
    println("Young BAM files found: ", length(bams2))
    exit(1)
end

# Append bam_dir path to each bam file
bams1 = [joinpath(bam_dir, bam) for bam in bams1]
bams2 = [joinpath(bam_dir, bam) for bam in bams2]

println("Old samples (", length(bams1), "):")
for bam in bams1
    println("  ", bam)
end

println("Young samples (", length(bams2), "):")
for bam in bams2
    println("  ", bam)
end

# Output paths
prefix = "cpeltdm"
outdir_final = "$(outdir)"

println("Starting CPEL-TDM analysis...")
println("Output directory: ", outdir_final)
println("Prefix: ", prefix)

# Run CPEL. If you have paired-end data, need to set pe=true
# If you have matched samples, need to set matched=true. 
# For matched samples, CPEL assumes that the two lists of bams
# are matched at the i-th index.
try
    cpel_tdm(bams1,bams2,bed,fa,outdir_final,prefix;pe=true,matched=false)
    println("CPEL-TDM analysis completed successfully!")
catch e
    println("ERROR in CPEL-TDM analysis: ", e)
    println("Stack trace:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
    exit(1)
end