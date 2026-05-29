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
        if filename.startswith("hg38_circrna_"):
            annotation_dict.setdefault("cscd", []).append(path)
        else:
            annotation_dict[filename.replace(".parquet", "")] = path
    return annotation_dict


def _download_required_files(
    s3,
    bucket: str,
    local_dir: str,
    other_files: List[str],
    cscd_files: List[str],
) -> Dict[str, object]:
    """Download required annotation tables and return file mapping."""
    annotation_dict: Dict[str, object] = {}

    for filename in other_files:
        object_key = f"db_tables/{filename}"
        local_path = os.path.join(local_dir, filename)
        try:
            s3.download_file(bucket, object_key, local_path)
            annotation_dict[filename.replace(".parquet", "")] = local_path
        except Exception as e:
            print(f"Warning: Failed to download {filename} -> {e}")

    for filename in cscd_files:
        object_key = f"CSCD_cleaned/{filename}"
        local_path = os.path.join(local_dir, filename)
        try:
            s3.download_file(bucket, object_key, local_path)
            annotation_dict.setdefault("cscd", []).append(local_path)
        except Exception as e:
            print(f"Warning: Failed to download {filename} -> {e}")

    return annotation_dict


def fetch_annotation_tables(lookup_results: Dict[str, Dict[str, pl.DataFrame]], tmp_dir_path: str = "tmp"):
    """
    Check if we need all CSCD annotations, or just a subset based on the lookup results.
    """
    cscd_chrs = _extract_cscd_chromosomes(lookup_results)
    cscd_files = [f"hg38_circrna_{chrom}.parquet" for chrom in cscd_chrs]
    other_files = [
        'arraystar.parquet',
        'circbank.parquet',
        'circbase.parquet',
        'circpedia.parquet',
        'circRNA_DB.parquet'
    ]

    required_files = other_files + cscd_files

    annotation_sums = os.path.join(os.getcwd(), "assets", "annotation_md5sum.csv")
    tmp_dir = os.path.join(os.getcwd(), tmp_dir_path)
    valid_paths = check_sums(
        tmp_dir=tmp_dir,
        md5sum_file=annotation_sums,
        dir_prefix="annotation_tables_",
        required_files=required_files,
    )
    if valid_paths:
        return _build_annotation_dict_from_paths(valid_paths)

    s3 = boto3.client(
        's3',
        region_name='eu-north-1',
        config=Config(signature_version=UNSIGNED)
    )

    bucket = 'digbyb'
    
    cwd_tmp = os.path.join(os.getcwd(), tmp_dir_path)
    os.makedirs(cwd_tmp, exist_ok=True)
    local_dir = tempfile.mkdtemp(prefix="annotation_tables_", dir=cwd_tmp)
    annotation_dict = _download_required_files(
        s3=s3,
        bucket=bucket,
        local_dir=local_dir,
        other_files=other_files,
        cscd_files=cscd_files,
    )

    cscd_count = len(annotation_dict.get("cscd", []))
    non_cscd_count = len([k for k in annotation_dict.keys() if k != "cscd"])
    dl_files = non_cscd_count + cscd_count
    print(f"Downloaded {dl_files} annotation files to {local_dir}")

    return annotation_dict


