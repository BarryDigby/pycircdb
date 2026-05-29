import os
import tempfile
import boto3
from botocore import UNSIGNED
from botocore.config import Config
from typing import List
from utils.md5sum_check import check_sums


def fetch_rbp_tables(tmp_dir_path: str = "tmp") -> List[str]:
    """
    Download (or return cached) RBP chromosome parquet files.
    """
    rbp_files = (
        [f"hg38_rbp_chr{chrom}.parquet" for chrom in range(1, 23)]
        + ["hg38_rbp_chrX.parquet", "hg38_rbp_chrY.parquet"]
    )

    # Check whether cached files in tmp/ already match the expected checksums.
    rbp_sums = os.path.join(os.getcwd(), "assets", "rbp_md5sum.csv")
    tmp_path = os.path.join(os.getcwd(), tmp_dir_path)
    valid_paths = check_sums(
        tmp_dir=tmp_path,
        md5sum_file=rbp_sums,
        dir_prefix="rbp_tables_",
        required_files=rbp_files,
    )
    if valid_paths:
        return valid_paths

    # Configure boto3 for anonymous access (--no-sign-request)
    s3 = boto3.client(
        's3', 
        region_name='eu-north-1',
        config=Config(signature_version=UNSIGNED)
    )
    
    bucket = 'digbyb'
    
    os.makedirs(tmp_path, exist_ok=True)
    local_dir = tempfile.mkdtemp(prefix="rbp_tables_", dir=tmp_path)
    downloaded: List[str] = []

    for filename in rbp_files:
        object_key = f"RBP/{filename}"
        local_path = os.path.join(local_dir, filename)

        try:
            s3.download_file(bucket, object_key, local_path)
            downloaded.append(local_path)
        except Exception as e:
            print(f"Warning: Failed to download {filename} -> {e}")

    print(f"Downloaded {len(downloaded)} RBP files to {local_dir}")
    return downloaded