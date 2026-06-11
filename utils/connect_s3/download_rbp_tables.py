import os
import boto3
import polars as pl
from botocore import UNSIGNED
from botocore.config import Config
from typing import List, Dict
from utils.md5sum_check import _load_expected_sums, _file_md5sum

def _extract_cscd_chromosomes(lookup_results: Dict[str, Dict[str, pl.DataFrame]]) -> List[str]:
    """Collect required CSCD chromosome names from lookup hits."""
    cscd_chrs = set()

    for _, db_dict in lookup_results.items():
        cscd_hits = db_dict.get("cscd")
        if cscd_hits is None or cscd_hits.shape == (0, 0):
            continue
        if "hg19" not in cscd_hits.columns or "hg38" not in cscd_hits.columns:
            continue

        cscd_hits = cscd_hits.with_columns(
            pl.col("hg19")
            .str.splitn(":", 2)
            .struct.rename_fields(["19chroms", "19rest"])
            .alias("hg19_fields")
        ).unnest("hg19_fields")

        cscd_hits = cscd_hits.with_columns(
            pl.col("hg38")
            .str.splitn(":", 2)
            .struct.rename_fields(["38chroms", "38rest"])
            .alias("hg38_fields")
        ).unnest("hg38_fields")

        cscd_chrs.update(cscd_hits["38chroms"].unique().to_list())
        cscd_chrs.update(cscd_hits["19chroms"].unique().to_list())

    cscd_chrs -= {None}
    return sorted(cscd_chrs)

def fetch_rbp_tables(lookup_results: Dict[str, Dict[str, pl.DataFrame]], tmp_dir_path: str = "tmp") -> List[str]:
    """
    Download (or return cached) RBP chromosome parquet files based on lookup results.
    """
    cscd_chrs = _extract_cscd_chromosomes(lookup_results)
    rbp_files = [f"hg38_rbp_{chrom}.parquet" for chrom in cscd_chrs]

    rbp_sums = os.path.join(os.getcwd(), "assets", "rbp_md5sum.csv")
    expected_sums = _load_expected_sums(rbp_sums)
    
    local_dir = os.path.join(os.getcwd(), tmp_dir_path, "rbp_tables")
    os.makedirs(local_dir, exist_ok=True)
    
    missing_files = []
    valid_paths = []
    
    for filename in rbp_files:
        file_path = os.path.join(local_dir, filename)
        if os.path.isfile(file_path) and _file_md5sum(file_path).lower() == expected_sums.get(filename):
            valid_paths.append(file_path)
        else:
            missing_files.append(filename)

    if missing_files:
        s3 = boto3.client(
            's3', 
            region_name='eu-north-1',
            config=Config(signature_version=UNSIGNED)
        )
        bucket = 'digbyb'
        
        print(f"Downloading {len(missing_files)} missing RBP files to {local_dir}")
        for filename in missing_files:
            object_key = f"RBP/{filename}"
            local_path = os.path.join(local_dir, filename)
            
            try:
                s3.download_file(bucket, object_key, local_path)
                valid_paths.append(local_path)
            except Exception as e:
                print(f"Warning: Failed to download {filename} -> {e}")
    else:
        print(f"Using cached files from {local_dir}; all MD5 checks passed.")

    return valid_paths