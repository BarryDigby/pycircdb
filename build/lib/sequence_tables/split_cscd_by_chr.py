import polars as pl
from pathlib import Path

# Path to the input parquet file
data_path = "cscd_sequence.parquet"
# Output filename template
output_template = "hg38_sequence_chr{chr}.parquet"

# Chromosomes to split by
chromosomes = [str(i) for i in range(1, 23)] + ["M", "X", "Y"]

# Read the input parquet file
print(f"Reading {data_path} ...")
df = pl.read_parquet(data_path)


# For each chromosome, filter and write to a new parquet file using polars string methods
for chr_name in chromosomes:
    prefix = f"chr{chr_name}"
    # Accepts chrN, chrN:, chrN_
    mask = df["circRNA"].str.starts_with(prefix)

    chr_df = df.filter(mask)
    if chr_df.height > 0:
        out_path = output_template.format(chr=chr_name)
        print(f"Writing {out_path} with {chr_df.height} rows ...")
        chr_df.write_parquet(out_path)
    else:
        print(f"No rows for chr{chr_name}")

print("Done.")
