from hamilton import driver
from hamilton.execution import executors
from hamilton.htypes import Collect, Parallelizable
from typing import Dict, Any, Union, List, Generator
from utils.connect_s3.download_annotation_tables import fetch_annotation_tables

import polars as pl
import utils.annotate.annotate_subdag as annotate_subdag


AnnotateBroadcast = Dict[str, Union[str, List[str], pl.DataFrame, None]]


def broadcast_annotation(
    config: Dict[str, Any],
    annotation_tables: Dict[str, Union[str, List[str]]],
    lookup_results: Dict[str, Dict[str, pl.DataFrame]],
) -> Parallelizable[AnnotateBroadcast]:
    """
    Broadcast lookup hits and annotation parquet paths to downstream nodes.

    Args:
        config: JSON format of the config file.
        lookup_results: Lookup step results for each sample and database.

    Yields:
        A dictionary with sample metadata, annotation path(s), and lookup hits.
    """
    for sample_name, lookup_hits in lookup_results.items():
        for db_name, pl_hits in lookup_hits.items():
            annotation_parquet_path = annotation_tables.get(db_name)
            if pl_hits.is_empty():
                continue

            yield {
                "sample_name": sample_name,
                "output_dir": config['global_parameters'].get("output_dir"),
                "db_name": db_name,
                "annotation_parquet_path": annotation_parquet_path,
                "lookup_hits": pl_hits,
            }


def instatiate_annotation_subdag(
    broadcast_annotation: AnnotateBroadcast
) -> None:
    """
    For each sample, pop off the db_name and convert it to config dict.
    Hamilton within Hamilton to extract params and cast as config for subdag.
    """
    db_name = broadcast_annotation.get("db_name")
    per_sample_annotation = dict(broadcast_annotation)
    print(f'starting {per_sample_annotation.get("sample_name")} -- {db_name} annotation')

    dr = (
        driver.Builder()
        .with_config({'db_name': db_name})
        .with_modules(annotate_subdag)
        .build()
    )

    dr.execute(
        ['write_to_output_dir'],
        inputs={'per_sample_annotation': per_sample_annotation}
    )['write_to_output_dir']


def close_annotation(
    instatiate_annotation_subdag: Collect[Any]
) -> None:
    """
    Closes the annotaton subdag.
    """
    print("Annotation subdag complete.")
