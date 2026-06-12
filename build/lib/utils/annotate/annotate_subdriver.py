import polars as pl
from hamilton import driver
from hamilton.execution import executors
from typing import Dict, Any, Union, List
from hamilton.htypes import Parallelizable, Collect

import utils.annotate.annotate_subdag as annotate_subdag

def broadcast_lookup(
    config: Dict[str, Any],
    annotation_tables: Dict[str, Union[str, List[str]]],
    lookup_results: Dict[str, Dict[str, pl.DataFrame]],
) -> Parallelizable[Dict[str, Union[str, List[str], pl.DataFrame, None]]]:
    """
    Broadcast lookup hits and annotation parquet paths to downstream nodes.

    Args:
        config: JSON format of the config file.
        annotation_tables: Mapping of annotation database name to parquet path(s).
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
                "output_dir": config.get("output_dir"),
                "db_name": db_name,
                "annotation_parquet_path": annotation_parquet_path,
                "lookup_hits": pl_hits,
            }


def instatiate_annotation_subdag(
    broadcast_lookup: Dict[str, Union[str, List[str], pl.DataFrame, None]]
) -> None:
    """
    For each sample, pop off the db_name and convert it to config dict.
    Hamilton within Hamilton to extract params and cast as config for subdag.
    """
    print(f"Received broadcast_lookup: {broadcast_lookup}")
    db_name = broadcast_lookup.get("db_name")
    per_sample_annotation = dict(broadcast_lookup)
    per_sample_annotation.pop("db_name", None)

    dr = (
        driver.Builder()
        .enable_dynamic_execution(allow_experimental_mode=True)
        .with_config({'db_name': db_name})
        .with_modules(annotate_subdag)
        .build()
    )

    result = dr.execute(
        ['write_to_output_dir'],
        inputs={'per_sample_annotation': per_sample_annotation}
    )['write_to_output_dir']

    return result


def return_collected_results(instatiate_annotation_subdag: Collect[Any]) -> Any:
    return None