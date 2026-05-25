from hamilton.htypes import Collect, Parallelizable
from typing import Dict, Any, Union, List
from pathlib import Path
import os
import polars as pl
import gzip

from utils.mirna.mirna_driver import broadcast_mirna


RBPBroadcast = Dict[str, Union[str, str, pl.DataFrame, List, None]]


def broadcast_rbp(
    config: Dict[str, Any],    
    lookup_dict: Dict[str, Dict[str, pl.DataFrame]],
    rbp_tables: List[str]
) -> Parallelizable[RBPBroadcast]:
    """
    Broadcast lookup hits and RBP table paths to downstream nodes.

    Args:
        config: JSON format of the config file.
        lookup_dict: Lookup step results for each sample and database.
        rbp_tables: List of RBP parquet paths.

    Yields:
        Per database hit, per sample, RBP annotation.
    """
    for sample_name, lookup_hits in lookup_dict.items():
        for db_name, pl_hits in lookup_hits.items():
            for chromosome_rbp_path in rbp_tables:

                foo =  {
                    "sample_name": sample_name,
                    "output_dir": config['global_parameters'].get("output_dir"),
                    "db_name": db_name,
                    "lookup_hits": pl_hits,
                    "rbp_table": chromosome_rbp_path,
                }

                print(foo)
                yield foo


def rbp_hits(broadcast_rbp: RBPBroadcast) -> None:
    """
    No need for a driver, per chromosome processing handled in parallel here.
    """
    output_dir = broadcast_rbp["output_dir"]
    sample_name = broadcast_rbp["sample_name"]
    lookup_pl = broadcast_rbp["lookup_hits"]
    rbp_table = broadcast_rbp["rbp_table"]
    chromosome = Path(rbp_table).stem.split("_")[-1]
        
    hg38_ids = lookup_pl["hg38"].str.split("|").list.first()

    query = (
        pl.scan_parquet(rbp_table)
        .filter(pl.col("circRNA").is_in(hg38_ids))
        )

    df = query.collect(streaming=True)

    print(df)
    
    if not df.is_empty():
        p = Path(output_dir)
        if p.is_absolute():
            output_path = p / f"{sample_name}/hg38_{chromosome}_rbp_hits.txt.gz"
        else:
            cwd_tmp = os.path.join(os.getcwd(), output_dir, sample_name, f"hg38_{chromosome}_rbp_hits.txt.gz")
            output_path = Path(cwd_tmp)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.write_csv(output_path, separator='\t', include_header=True, compression='gzip')





def close_rbp(
    rbp_hits: Collect[Any]
) -> None:
    """
    Closes the RBP subdag.
    """
    print("RBP subdag complete.")