import os
import tempfile
import boto3
import polars as pl
from typing import Dict, List
from botocore import UNSIGNED
from botocore.config import Config

from utils.md5sum_check import check_sums

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

def _build_annotation_dict_from_paths(valid_paths: List[str]) -> Dict[str, object]:
    """Build annotation mapping from validated local file paths."""
    annotation_dict: Dict[str, object] = {}
    for path in valid_paths:
        filename = os.path.basename(path)
        if filename.startswith("hg38_sequence_"):
            annotation_dict.setdefault("cscd", []).append(path)
        else:
            annotation_dict[filename.replace(".parquet", "")] = path
    return annotation_dict


def _download_required_files(
    s3,
    bucket: str,
    local_dir: str,
    other_files: List[str],
    cscd_chroms: List[str]
) -> Dict[str, object]:
    """Download required sequence parquet files from S3, including CSCD chromosome-specific tables."""
    sequence_dict: Dict[str, object] = {}

    for filename in other_files:
        object_key = f"sequence_tables/{filename}"
        local_path = os.path.join(local_dir, filename)
        
        try:
            s3.download_file(bucket, object_key, local_path)
            sequence_dict[filename.replace("_sequence.parquet", "")] = local_path
        except Exception as e:
            print(f"Warning: Failed to download {filename} -> {e}")

    for chrom in cscd_chroms:
        cscd_filename = f"hg38_sequence_{chrom}.parquet"
        object_key = f"sequence_tables/{cscd_filename}"
        local_path = os.path.join(local_dir, cscd_filename)

        try:
            s3.download_file(bucket, object_key, local_path)
            sequence_dict.setdefault("cscd", []).append(local_path)
        except Exception as e:
            print(f"Warning: Failed to download {cscd_filename} -> {e}")
    
    return sequence_dict

def fetch_sequence_tables(lookup_results: Dict[str, Dict[str, pl.DataFrame]]):
    """
    1.5GB locally
    """
    cscd_chrs = _extract_cscd_chromosomes(lookup_results)
    cscd_files = [f"hg38_sequence_{chrom}.parquet" for chrom in cscd_chrs]
    other_files = [
        'arraystar_sequence.parquet',
        'circatlas_sequence.parquet',
        'circbank_sequence.parquet',
        'circbase_sequence.parquet',
        'circpedia_sequence.parquet',
        'circRNA_DB_sequence.parquet'
    ]

    required_files = other_files + cscd_files
    sequence_sums = os.path.join(os.getcwd(), "assets", "sequence_md5sum.csv")
    tmp_dir = os.path.join(os.getcwd(), "tmp")
    valid_paths = check_sums(
        tmp_dir=tmp_dir, 
        md5sum_file=sequence_sums, 
        dir_prefix="sequence_tables_",
        required_files=required_files
    )

    if valid_paths:
        return _build_annotation_dict_from_paths(valid_paths)    
    
    s3 = boto3.client(
        's3', 
        region_name='eu-north-1',
        config=Config(signature_version=UNSIGNED)
    )

    bucket = 'digbyb'

    cwd_tmp = os.path.join(os.getcwd(), "tmp")
    os.makedirs(cwd_tmp, exist_ok=True)
    local_dir = tempfile.mkdtemp(prefix="sequence_tables_", dir=cwd_tmp)
    sequence_dict = _download_required_files(s3, bucket, local_dir, other_files, cscd_chrs)
    
    cscd_count = len(sequence_dict.get("cscd", []))
    non_cscd_count = len(sequence_dict) - (1 if "cscd" in sequence_dict else 0)
    dl_files = non_cscd_count + cscd_count
    print(f"Downloaded {dl_files} sequence tables to {local_dir}")

    return sequence_dict