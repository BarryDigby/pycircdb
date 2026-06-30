import os
import boto3
import polars as pl
from botocore import UNSIGNED
from botocore.config import Config
from typing import List, Dict
from rich.console import Console
from rich.text import Text

from utils.md5sum_check import _load_expected_sums, _file_md5sum

console = Console(stderr=True, highlight=False)

def _extract_required_chromosomes(lookup_results: Dict[str, Dict[str, pl.DataFrame]]) -> List[str]:
    """Collect required hg38 chromosome names from hits across all lookup tables.

    miRNA hits are matched against the hg38 IDs of every database hit (see
    broadcast_mirna), so the chromosome tables we download must be derived from
    all lookup tables, not just CSCD.
    """
    chroms = set()

    for db_dict in lookup_results.values():
        for pl_hits in db_dict.values():
            if pl_hits is None or pl_hits.is_empty():
                continue
            if "hg38" not in pl_hits.columns:
                continue

            chrom_series = (
                pl_hits["hg38"]
                .str.split("|").list.first()   # drop optional strand suffix
                .str.split(":").list.first()   # chromosome is the leading field
            )
            chroms.update(chrom_series.unique().to_list())

    chroms -= {None}
    return sorted(chroms)

def fetch_mirna_tables(lookup_results: Dict[str, Dict[str, pl.DataFrame]], tmp_dir_path: str = "tmp", verbose: int = 1) -> List[str]:
    """
    Download (or return cached) miRNA chromosome parquet files based on lookup results.
    """
    required_chrs = _extract_required_chromosomes(lookup_results)
    mirna_files = [f"hg38_mirna_{chrom}.parquet" for chrom in required_chrs]

    from pathlib import Path
    assets_dir = Path(__file__).resolve().parent.parent.parent / "assets"
    mirna_sums = str(assets_dir / "mirna_md5sum.csv")
    expected_sums = _load_expected_sums(mirna_sums)
    
    local_dir = os.path.join(os.getcwd(), tmp_dir_path, "mirna_tables")
    os.makedirs(local_dir, exist_ok=True)
    
    if verbose >= 1:
        console.print(Text(f"Checking for miRNA tables in {local_dir}...", style="bold white"))
        
    missing_files = []
    valid_paths = []
    
    for filename in mirna_files:
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
        
        if verbose >= 1:
            console.print(Text(f"Downloading {len(missing_files)} missing miRNA files to {local_dir}...", style="bold yellow"))
            
        for filename in missing_files:
            if verbose >= 2:
                console.print(f"  Downloading {filename}", style="cyan")
            object_key = f"miRNA/{filename}"
            local_path = os.path.join(local_dir, filename)
            
            try:
                s3.download_file(bucket, object_key, local_path)
                valid_paths.append(local_path)
            except Exception as e:
                console.print(f"Warning: Failed to download {filename} -> {e}", style="bold red")
                
        if verbose >= 1:
            console.print(Text("Successfully downloaded missing miRNA tables.", style="bold green"))
    else:
        if verbose >= 1:
            console.print(Text(f"Using cached files from {local_dir}; all MD5 checks passed.", style="bold green"))

    return valid_paths