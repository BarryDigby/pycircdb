import polars as pl
from typing import Dict, Any, Tuple
from polars import col
from hamilton.htypes import Parallelizable, Collect


def broadcast_config(config: Dict[str, Any], lookup_tables: Dict[str, str]) -> Parallelizable[Dict[str, Any]]:
    """
    Load Parquet Files once, then broadcast with config params to next node.
    
    Args:
        config: JSON format of the config file.
        lookup_tables: Dictionary mapping database names to their respective Parquet file paths.
    
    Yields:
        A dictionary containing sample name, user input file path, parameters and the loaded
        lookup table for each database, which will be processed in parallel in the next node.
    """
    for db_name, lookup_path in lookup_tables.items():
        lookup_pl = pl.read_parquet(lookup_path)
        for sample_name, sample_info in config.get("samples", {}).items():
            combined_input = {
                "sample_name": sample_name,
                "file_path": sample_info.get("input"),
                "reference": sample_info.get("reference"),
                "zero_based": sample_info.get("zero_based"),
                "lookup_pl": lookup_pl,
                "db_name": db_name
            }
            yield combined_input

def database_lookup(broadcast_config: Dict[str, Any]) -> Dict[str, Dict[str, pl.DataFrame]]:
    """
    Perform simple filter to return {sample_name: {db_name: hits}}
    
    Returns:
        A dictionary for each sample containing a dictionary with the database name and the corresponding hits DataFrame.
    """
    sample_name = broadcast_config.get("sample_name")
    reference = broadcast_config.get("reference")
    lookup_pl = broadcast_config.get("lookup_pl")
    db_name = broadcast_config.get("db_name")
    zero_based = broadcast_config.get("zero_based", True)
    
    input_pl = pl.read_csv(broadcast_config.get("file_path"), has_header=False, new_columns=['ID'])
    
    if not zero_based:
        input_pl = input_pl.with_columns(
            pl.when(col('ID').str.contains(r'^.*?:\d+-\d+.*$'))
            .then(
                pl.concat_str([
                    col('ID').str.extract(r'^(.*?:)(\d+)(-.*)$', 1),
                    (col('ID').str.extract(r'^(.*?:)(\d+)(-.*)$', 2).cast(pl.Int64) - 1).cast(pl.Utf8),
                    col('ID').str.extract(r'^(.*?:)(\d+)(-.*)$', 3)
                ])
            )
            .otherwise(col('ID'))
            .alias('ID')
        )

    hits = lookup_pl.filter(col(f"{reference}").is_in(input_pl['ID']))

    #print(f"Hits for sample {sample_name} in database {db_name}:")
    #print(hits)
    
    return {sample_name: {db_name: hits}}


def return_collected_results(database_lookup: Collect[Dict[str, Dict[str, pl.DataFrame]]]) -> Dict[str, Dict[str, pl.DataFrame]]:
    """
    Collate sample results into a single dictionary of {sample_name: {db_name: DataFrame}}"
    
    Returns:
        A reduced dictionary where each sample name maps to a dictionary of database names and their corresponding hits DataFrames.
    """
    merged_results = {}
    for sample_dict in database_lookup:
        for sample_name, db_dict in sample_dict.items():
            if sample_name not in merged_results:
                merged_results[sample_name] = {}
            merged_results[sample_name].update(db_dict)
    return merged_results
