import pandas as pd

# Load CSV
df = pd.read_csv("/dcs07/afeinber/data/personal/jzhan/CPEL2/DMR_old-vs-young_group2_corrected.csv", index_col=0)

# Create BED DataFrame
bed = pd.DataFrame({
    "chrom": df["chr"],
    "start": df["start"].astype(int),
    "end": df["end"].astype(int),
    "name": ["DMR" + str(idx) for idx in df.index],
    "score": df["areaStat"].round(3),  # or use "0" if no score
    "strand": "."  # No strand info
})

# Save as BED (tab-separated, no header)
bed.to_csv("output.bed", sep="\t", header=False, index=False)