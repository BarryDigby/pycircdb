from hamilton.htypes import Collect, Parallelizable
from typing import Dict, Any, Union, List
from pathlib import Path
from rich.console import Console
import os
import polars as pl

console = Console(stderr=True, highlight=False)

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
        # Combine all hg38 IDs for this sample across all databases
        hg38_series = []
        # Map strand-stripped hg38 key -> stranded hg38/hg19 coordinates,
        # so both can be restored on the output after the strand-less match.
        strand_frames = []
        for db_name, pl_hits in lookup_hits.items():
            if pl_hits.is_empty() or "hg38" not in pl_hits.columns:
                continue
            hg38_series.append(pl_hits["hg38"].str.split("|").list.first())
            hg19_col = pl.col("hg19") if "hg19" in pl_hits.columns else pl.lit(None, dtype=pl.Utf8).alias("hg19")
            strand_frames.append(pl_hits.select(
                pl.col("hg38").str.split("|").list.first().alias("_hg38_key"),
                pl.col("hg38"),
                hg19_col,
            ))

        if not hg38_series:
            continue

        unique_hg38_ids = pl.concat(hg38_series).unique().to_list()

        # Prefer stranded coordinates when de-duplicating the strand-less key.
        strand_map = (
            pl.concat(strand_frames)
            .sort("hg38")
            .unique(subset="_hg38_key", keep="last")
        )

        for chromosome_rbp_path in rbp_tables:
            yield {
                "sample_name": sample_name,
                "output_dir": config['global_parameters'].get("output_dir"),
                "unique_hg38_ids": unique_hg38_ids,
                "strand_map": strand_map,
                "rbp_table": chromosome_rbp_path,
                "verbose": config.get("verbose", 1)
            }



def rbp_hits(broadcast_rbp: RBPBroadcast) -> None:
    """
    No need for a driver, per chromosome processing handled in parallel here.
    """
    output_dir = broadcast_rbp["output_dir"]
    sample_name = broadcast_rbp["sample_name"]
    unique_hg38_ids = broadcast_rbp["unique_hg38_ids"]
    strand_map = broadcast_rbp.get("strand_map")
    rbp_table = broadcast_rbp["rbp_table"]
    verbose = broadcast_rbp.get("verbose", 1)
    chromosome = Path(rbp_table).stem.split("_")[-1]
        
    if verbose >= 2:
        console.print(f"  Starting RBP extraction for {sample_name}: [cyan]processing {chromosome}[/cyan]")

    query = (
        pl.scan_parquet(rbp_table)
        .filter(pl.col("circRNA").is_in(unique_hg38_ids))
        )

    df = query.collect(engine='streaming')
    
    if not df.is_empty():
        # The RBP table's circRNA column is the strand-stripped hg38 key.
        # Join back to restore stranded hg38/hg19 coordinates onto the output.
        if strand_map is not None and strand_map.height > 0:
            df = (
                df.rename({"circRNA": "_hg38_key"})
                .join(strand_map, on="_hg38_key", how="left")
                # Fall back to the strand-less key if hg38 is unmapped.
                .with_columns(pl.col("hg38").fill_null(pl.col("_hg38_key")))
            )
        else:
            # No strand info available: fall back to the strand-less key for hg38.
            df = df.rename({"circRNA": "hg38"}).with_columns(pl.lit(None, dtype=pl.Utf8).alias("hg19"))

        ordered = ["hg19", "hg38"] + [c for c in df.columns if c not in ("hg19", "hg38", "_hg38_key")]
        df = df.select(ordered)

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
    console.print("RBP subdag complete.", style="bold green")