from hamilton.htypes import Collect, Parallelizable
from typing import Dict, Any, Union, List
from pathlib import Path
import os
import polars as pl
import gzip


miRNABroadcast = Dict[str, Union[str, str, pl.DataFrame, List, None]]


def broadcast_mirna(
    config: Dict[str, Any],    
    lookup_dict: Dict[str, Dict[str, pl.DataFrame]],
    mirna_tables: List[str]
) -> Parallelizable[miRNABroadcast]:
    """
    Broadcast lookup hits and miRNA table paths to downstream nodes.

    Args:
        config: JSON format of the config file.
        lookup_dict: Lookup step results for each sample and database.
        mirna_tables: List of miRNA parquet paths.

    Yields:
        Per database hit, per sample, miRNA annotation.
    """
    for sample_name, lookup_hits in lookup_dict.items():
        for db_name, pl_hits in lookup_hits.items():
            for chromosome_mirna_path in mirna_tables:

                foo =  {
                    "sample_name": sample_name,
                    "output_dir": config['global_parameters'].get("output_dir"),
                    "db_name": db_name,
                    "lookup_hits": pl_hits,
                    "mirna_table": chromosome_mirna_path,
                }

                print(foo)
                yield foo


def mirna_hits(broadcast_mirna: miRNABroadcast) -> None:
    """
    No need for a driver, per chromosome processing handled in parallel here.
    """
    output_dir = broadcast_mirna["output_dir"]
    sample_name = broadcast_mirna["sample_name"]
    lookup_pl = broadcast_mirna["lookup_hits"]
    mirna_table = broadcast_mirna["mirna_table"]
    chromosome = Path(mirna_table).stem.split("_")[-1]
        
    hg38_ids = lookup_pl["hg38"].str.split("|").list.first()

    query = (
        pl.scan_parquet(mirna_table)
        .filter(pl.col("circRNA").is_in(hg38_ids))
        )

    df = query.collect(streaming=True)

    print(df)
    
    if not df.is_empty():
        p = Path(output_dir)
        if p.is_absolute():
            output_path = p / f"{sample_name}/hg38_{chromosome}_mirna_hits.txt.gz"
        else:
            cwd_tmp = os.path.join(os.getcwd(), output_dir, sample_name, f"hg38_{chromosome}_mirna_hits.txt.gz")
            output_path = Path(cwd_tmp)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.write_csv(output_path, separator='\t', include_header=True, compression='gzip')





def close_mirna(
    mirna_hits: Collect[Any]
) -> None:
    """
    Closes the miRNA subdag.
    """
    print("miRNA subdag complete.")