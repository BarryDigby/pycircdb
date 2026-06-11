import polars as pl
from polars import col
from hamilton.function_modifiers import config
from typing import Dict, Any, Union, List
from pathlib import Path
import os


PerSampleSequence = Dict[str, Union[str, str, pl.DataFrame, List, None]]


def parquet_to_fasta(parquet_table: pl.DataFrame, output_fasta_path: str) -> None:
    """
    Helper function, not within Hamilton scope. 
    Convert a parquet file with two columns to FASTA format.
    """
    with open(output_fasta_path, 'w') as f:
        for row in parquet_table.iter_rows():
            identifier, sequence = row[0], row[1]
            f.write(f">{identifier}\n{sequence}\n")


@config.when(db_name='arraystar')
def sequences__arraystar(
    per_sample_sequence: PerSampleSequence
) -> Dict[str, Dict[str, pl.DataFrame]]:
    """
    arraystar, sequence
    """
    lookup_pl = per_sample_sequence.get("lookup_hits")

    query = (
        pl.scan_parquet(per_sample_sequence.get("sequence_tables"))
        .filter(col("arraystar").is_in(lookup_pl['arraystar']))
    )

    df = query.collect(streaming=True)

    return {per_sample_sequence.get("sample_name"): {'arraystar': df}}


@config.when(db_name='circatlas')
def sequences__circatlas(
    per_sample_sequence: PerSampleSequence
) -> Dict[str, Dict[str, pl.DataFrame]]:
    """
    circatlas, sequence
    """
    lookup_pl = per_sample_sequence.get("lookup_hits")

    query = (
        pl.scan_parquet(per_sample_sequence.get("sequence_tables"))
        .filter(col("circatlas").is_in(lookup_pl['circatlas']))
    )

    df = query.collect(streaming=True)

    return {per_sample_sequence.get("sample_name"): {'circatlas': df}}


@config.when(db_name='circbank')
def sequences__circbank(
    per_sample_sequence: PerSampleSequence
) -> Dict[str, Dict[str, pl.DataFrame]]:
    """
    circbank, sequence
    """
    lookup_pl = per_sample_sequence.get("lookup_hits")

    query = (
        pl.scan_parquet(per_sample_sequence.get("sequence_tables"))
        .filter(col("circbank").is_in(lookup_pl['circbank']))
    )

    df = query.collect(streaming=True)

    return {per_sample_sequence.get("sample_name"): {'circbank': df}}


@config.when(db_name='circbase')
def sequences__circbase(
    per_sample_sequence: PerSampleSequence
) -> Dict[str, Dict[str, pl.DataFrame]]:
    """
    circbase, sequence
    """
    lookup_pl = per_sample_sequence.get("lookup_hits")

    query = (
        pl.scan_parquet(per_sample_sequence.get("sequence_tables"))
        .filter(col("circbase").is_in(lookup_pl['circbase']))
    )

    df = query.collect(streaming=True)

    return {per_sample_sequence.get("sample_name"): {'circbase': df}}


@config.when(db_name='circpedia')
def sequences__circpedia(
    per_sample_sequence: PerSampleSequence
) -> Dict[str, Dict[str, pl.DataFrame]]:
    """
    circpedia, sequence
    """
    lookup_pl = per_sample_sequence.get("lookup_hits")

    query = (
        pl.scan_parquet(per_sample_sequence.get("sequence_tables"))
        .filter(col("circpedia").is_in(lookup_pl['circpedia']))
    )

    df = query.collect(streaming=True)

    return {per_sample_sequence.get("sample_name"): {'circpedia': df}}


@config.when(db_name='circRNA_DB')
def sequences__circRNA_DB(
    per_sample_sequence: PerSampleSequence
) -> Dict[str, Dict[str, pl.DataFrame]]:
    """
    circRNA_DB, sequence
    """
    lookup_pl = per_sample_sequence.get("lookup_hits")

    query = (
        pl.scan_parquet(per_sample_sequence.get("sequence_tables"))
        .filter(col("circRNA_DB").is_in(lookup_pl['circRNA_DB']))
    )

    df = query.collect(streaming=True)

    return {per_sample_sequence.get("sample_name"): {'circRNA_DB': df}}


@config.when(db_name='cscd')
def sequences__cscd(
    per_sample_sequence: PerSampleSequence
) -> Dict[str, Dict[str, pl.DataFrame]]:
    """
    circRNA, sequence
    """
    lookup_pl = per_sample_sequence.get("lookup_hits")
    sequence_paths = per_sample_sequence.get("sequence_tables")
    
    queries = []
    for sequence_path in sequence_paths:
        query = (
            pl.scan_parquet(sequence_path)
            .filter(col("circRNA").is_in(lookup_pl['hg38']))
        )
        queries.append(query)
        
    df = pl.concat(queries, how='vertical').collect(streaming=True)

    return {per_sample_sequence.get("sample_name"): {'cscd': df}}


def write_to_output_dir(
    sequences: Dict[str, Dict[str, pl.DataFrame]],
    per_sample_sequence: PerSampleSequence
) -> None:
    """
    Writes sequence hits to the output directory.
    """
    output_dir = per_sample_sequence.get("output_dir")
    sample_name = per_sample_sequence.get("sample_name")
    db_name = per_sample_sequence.get("db_name")

    sequence_hits = sequences.get(sample_name, {}).get(db_name)

    p = Path(output_dir).expanduser()
    if p.is_absolute():
        output_path = p / sample_name
    else:
        cwd_tmp = os.path.join(os.getcwd(), output_dir, sample_name)
        output_path = Path(cwd_tmp)

    output_path.mkdir(parents=True, exist_ok=True)
    output_file = output_path / f"{db_name}.fasta"
    
    parquet_to_fasta(sequence_hits, output_file)