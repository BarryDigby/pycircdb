from hamilton.htypes import Collect, Parallelizable
from typing import Dict, Any, Union, List
from pathlib import Path
from rich.console import Console
import os
import polars as pl

console = Console(stderr=True, highlight=False)


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
        # Combine all hg38 IDs for this sample across all databases
        hg38_series = []
        # Map strand-stripped hg38 key -> full (stranded) hg38 coordinate,
        # so strand can be restored on the output after the strand-less match.
        strand_frames = []
        for db_name, pl_hits in lookup_hits.items():
            if pl_hits.is_empty() or "hg38" not in pl_hits.columns:
                continue
            hg38_series.append(pl_hits["hg38"].str.split("|").list.first())
            strand_frames.append(pl_hits.select(
                pl.col("hg38").str.split("|").list.first().alias("_hg38_key"),
                pl.col("hg38").alias("circRNA"),
            ))

        if not hg38_series:
            continue

        unique_hg38_ids = pl.concat(hg38_series).unique().to_list()

        # Prefer stranded coordinates when de-duplicating the strand-less key.
        strand_map = (
            pl.concat(strand_frames)
            .sort("circRNA")
            .unique(subset="_hg38_key", keep="last")
        )

        for chromosome_mirna_path in mirna_tables:
            yield {
                "sample_name": sample_name,
                "output_dir": config['global_parameters'].get("output_dir"),
                "unique_hg38_ids": unique_hg38_ids,
                "strand_map": strand_map,
                "mirna_table": chromosome_mirna_path,
                "mirna_algorithms": config.get("mirna_algorithms"),
                "verbose": config.get("verbose", 1)
            }


def mirna_hits(broadcast_mirna: miRNABroadcast) -> None:
    """
    No need for a driver, per chromosome processing handled in parallel here.
    """
    output_dir = broadcast_mirna["output_dir"]
    sample_name = broadcast_mirna["sample_name"]
    unique_hg38_ids = broadcast_mirna["unique_hg38_ids"]
    strand_map = broadcast_mirna.get("strand_map")
    mirna_table = broadcast_mirna["mirna_table"]
    mirna_algorithms = broadcast_mirna.get("mirna_algorithms")
    verbose = broadcast_mirna.get("verbose", 1)
    chromosome = Path(mirna_table).stem.split("_")[-1]
    
    if verbose >= 2:
        console.print(f"  Starting miRNA extraction for {sample_name}: [cyan]processing {chromosome}[/cyan]")
        
    query = (
        pl.scan_parquet(mirna_table)
        .filter(pl.col("circRNA").is_in(unique_hg38_ids))
    )

    if mirna_algorithms:
        # Use simple str.contains across the algorithms provided, case insensitive
        pattern = f"(?i){'|'.join(mirna_algorithms)}"
        query = query.filter(pl.col("Algorithm").str.contains(pattern))

    df = query.collect(engine="streaming")
    
    if not df.is_empty():
        # The miRNA table's circRNA column is the strand-stripped hg38 key.
        # Join back to restore the strand onto the circRNA coordinate.
        if strand_map is not None and strand_map.height > 0:
            df = (
                df.rename({"circRNA": "_hg38_key"})
                .join(strand_map, on="_hg38_key", how="left")
                # Fall back to the strand-less key if a coordinate is unmapped.
                .with_columns(pl.col("circRNA").fill_null(pl.col("_hg38_key")))
            )
            ordered = ["circRNA"] + [c for c in df.columns if c not in ("circRNA", "_hg38_key")]
            df = df.select(ordered)

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
    console.print("miRNA subdag complete.", style="bold green")