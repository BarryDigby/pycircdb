#!/usr/bin/env python3
"""
Extract sequences from hg38_circrna parquet files and save separately.
"""
import polars as pl
from pathlib import Path
from glob import glob

def main():
    # Get current working directory
    cwd = Path.cwd()
    
    # Create list of all matching parquet files
    chromosomes = list(range(1, 23)) + ['X', 'Y', 'M']
    parquet_files = []
    
    for chrom in chromosomes:
        pattern = f"hg38_circrna_chr{chrom}.parquet"
        matching = glob(str(cwd / pattern))
        parquet_files.extend(matching)
    
    parquet_files.sort()
    
    if not parquet_files:
        print("No matching parquet files found.")
        return
    
    print(f"Found {len(parquet_files)} parquet files")
    
    # List to collect circRNA + sequence dataframes
    sequence_dfs = []
    
    # Process each file
    for filepath in parquet_files:
        print(f"Processing {Path(filepath).name}...")
        
        # Read parquet file
        df = pl.read_parquet(filepath)
        
        # Create new dataframe with only circRNA and sequence
        if "circRNA" in df.columns and "sequence" in df.columns:
            seq_df = df.select(["circRNA", "sequence"])
            sequence_dfs.append(seq_df)
            
            # Drop sequence from original dataframe and save back
            df_no_seq = df.drop("sequence")
            df_no_seq.write_parquet(filepath)
            print(f"  - Saved {Path(filepath).name} without sequence column")
        else:
            missing = []
            if "circRNA" not in df.columns:
                missing.append("circRNA")
            if "sequence" not in df.columns:
                missing.append("sequence")
            print(f"  WARNING: Missing columns {missing}, skipping...")
    
    # Combine all sequence dataframes into master
    if sequence_dfs:
        master_df = pl.concat(sequence_dfs)
        output_path = cwd / "cscd_sequence.parquet"
        master_df.write_parquet(output_path)
        print(f"\nMaster sequence file saved to {output_path}")
        print(f"Total rows: {len(master_df)}")
    else:
        print("No sequence dataframes were collected.")

if __name__ == "__main__":
    main()
