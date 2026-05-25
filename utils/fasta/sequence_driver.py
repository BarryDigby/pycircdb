import polars as pl
from hamilton import driver
from hamilton.execution import executors
from hamilton.htypes import Collect, Parallelizable
from typing import Dict, Any, Union, List
from utils.connect_s3.download_sequence_tables import fetch_sequence_tables
import utils.fasta.sequence_subdag as sequence_subdag

SequenceBroadcast = Dict[str, Union[str, str, pl.DataFrame, List, None]]

def broadcast_sequence(
    config: Dict[str, Any],
    lookup_dict: Dict[str, Dict[str, pl.DataFrame]],
    sequence_tables: Dict[str, str]
) -> Parallelizable[SequenceBroadcast]:
    """
    Broadcast lookup hits and sequence parquet paths to downstream nodes.

    Args:
        config: JSON format of the config file.
        lookup_dict: Lookup step results for each sample and database.
        sequence_tables: Dictionary of sequence table paths keyed by database name.
    Yields:
        A dictionary with sample metadata, sequence path(s), and lookup hits.
    """
    for sample_name, lookup_hits in lookup_dict.items():
        for db_name, pl_hits in lookup_hits.items():
            if pl_hits.is_empty():
                continue

            yield {
                "sample_name": sample_name,
                "output_dir": config['global_parameters'].get("output_dir"),
                "db_name": db_name,
                "lookup_hits": pl_hits,
                "sequence_tables": sequence_tables.get(db_name),
            }

def instantiate_sequence_subdag(broadcast_sequence: SequenceBroadcast) -> None:
    """
    Builds a Hamilton driver, passes sequence table paths and config to driver.
    """
    db_name = broadcast_sequence.get("db_name")
    per_sample_sequence = dict(broadcast_sequence)

    dr = (
        driver.Builder()
        .with_config({'db_name': db_name})
        .with_modules(sequence_subdag)
        .build()
    )

    dr.execute(
        ['write_to_output_dir'], 
        inputs={'per_sample_sequence': per_sample_sequence}
    )['write_to_output_dir']


def close_sequence(
    instantiate_sequence_subdag: Collect[Any]
) -> None:
    """
    Closes the sequence subdag.
    """
    print("Sequence subdag complete.")
