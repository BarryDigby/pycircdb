import os
import tempfile
import boto3
from botocore import UNSIGNED
from botocore.config import Config

from utils.md5sum_check import check_sums


def fetch_lookup_tables(tmp_dir_path: str = "tmp") -> Dict[str, str]:
    """
    80MB storage used locally.
    """
    known_files = [
        "arraystar_lookup.parquet",
        "circatlas_lookup.parquet",
        "circbank_lookup.parquet",
        "circbase_lookup.parquet",
        "circpedia_lookup.parquet",
        "circRNA_DB_lookup.parquet",
        "cscd_lookup.parquet"
    ]
    
    # Check whether cached files in tmp/ already match the expected checksums.
    lookup_sums = os.path.join(os.getcwd(), "assets", "lookup_md5sum.csv")
    tmp_dir = os.path.join(os.getcwd(), tmp_dir_path)
    valid_paths = check_sums(tmp_dir=tmp_dir, md5sum_file=lookup_sums, dir_prefix="lookup_tables_", required_files=known_files)
    if valid_paths:
        lookup_dict = {
            os.path.basename(path).replace("_lookup.parquet", ""): path
            for path in valid_paths
        }
        return lookup_dict

    # Configure boto3 for anonymous access (--no-sign-request)
    s3 = boto3.client(
        's3', 
        region_name='eu-north-1',
        config=Config(signature_version=UNSIGNED)
    )
    
    bucket = 'digbyb'
    
    # Create a temporary directory to store downloaded files
    cwd_tmp = os.path.join(os.getcwd(), tmp_dir_path)
    os.makedirs(cwd_tmp, exist_ok=True)
    local_dir = tempfile.mkdtemp(prefix="lookup_tables_", dir=cwd_tmp)
    lookup_dict = {}
    
    # Download explicitly requested files
    for filename in known_files:
        object_key = f"lookup_tables/{filename}"
        local_path = os.path.join(local_dir, filename)
        
        try:
            s3.download_file(bucket, object_key, local_path)
            key = filename.replace("_lookup.parquet", "")
            lookup_dict[key] = local_path
        except Exception as e:
            print(f"Warning: Failed to download {filename} -> {e}")

    print(f"Downloaded {len(lookup_dict)} lookup files to {local_dir}")
    
    # If the user doesn't want to keep it, you can delete it later using:
    # shutil.rmtree(local_dir)
    
    return lookup_dict