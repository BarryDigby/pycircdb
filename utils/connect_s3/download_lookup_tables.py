import os
import boto3
from typing import Dict
from botocore import UNSIGNED
from botocore.config import Config

from utils.md5sum_check import _load_expected_sums, _file_md5sum


def fetch_lookup_tables(tmp_dir_path: str = "tmp") -> Dict[str, str]:
    """
    85MB storage used locally.
    """
    known_files = [
        "arraystar_lookup.parquet",
        "circatlas_lookup.parquet",
        "circbank_lookup.parquet",
        "circbase_lookup.parquet",
        "circpedia_lookup.parquet",
        "circRNA_DB_lookup.parquet",
        "cscd_lookup.parquet",
        "exorbase_lookup.parquet"
    ]
    
    lookup_sums = os.path.join(os.getcwd(), "assets", "lookup_md5sum.csv")
    expected_sums = _load_expected_sums(lookup_sums)
    
    local_dir = os.path.join(os.getcwd(), tmp_dir_path, "lookup_tables")
    os.makedirs(local_dir, exist_ok=True)
    
    missing_files = []
    valid_paths = []
    
    for filename in known_files:
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
        
        print(f"Downloading {len(missing_files)} missing lookup files to {local_dir}")
        for filename in missing_files:
            object_key = f"lookup_tables/{filename}"
            local_path = os.path.join(local_dir, filename)
            
            try:
                s3.download_file(bucket, object_key, local_path)
                valid_paths.append(local_path)
            except Exception as e:
                print(f"Warning: Failed to download {filename} -> {e}")
    else:
        print(f"Using cached files from {local_dir}; all MD5 checks passed.")

    lookup_dict = {
        os.path.basename(path).replace("_lookup.parquet", ""): path
        for path in valid_paths
    }
    
    return lookup_dict