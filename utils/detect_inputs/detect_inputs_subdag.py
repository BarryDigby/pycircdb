import os
import re
import polars as pl
from pathlib import Path
from typing import Dict, Any, Tuple
from polars import col
from hamilton.htypes import Parallelizable, Collect


# Matches 'chr:start-end' with an optional '|+'/'|-' strand suffix.
_COORD_RE = re.compile(r'^(.*?):(\d+)-(\d+)(\|[+-])?$')


def _expand_input_keys(ids) -> Tuple[set, set]:
    """Expand input circRNA coordinates

    Idea here is to let users pass stranded & nonstranded coordinates.
    Next, we generate three coordinates for each, +/- 1 and the original.
    e.g:
    has strand (stranded_keys)
    chr4:1228198-1241519|- >>> chr4:1228197-1241519|-, chr4:1228198-1241519|-, chr4:1228199-1241519|-
    no strand (posonly_keys):
    chr4:1228198-1241519 >>> chr4:1228197-1241519, chr4:1228198-1241519, chr4:1228199-1241519

    Inputs that are not coordinates are ignored, as matching is performed
    against the hg19/hg38 coordinate column.

    Returns:
        (stranded_keys, posonly_keys) where stranded_keys are matched verbatim
        against the reference column and posonly_keys against a strand-stripped
        reference.
    """
    stranded_keys = set()
    posonly_keys = set()

    for raw in ids:
        if raw is None:
            continue
        match = _COORD_RE.match(raw)
        if match is None:
            continue
        chrom, start, end, strand = match.group(1), int(match.group(2)), match.group(3), match.group(4)
        for shifted_start in (start, start - 1, start + 1):
            if shifted_start < 0:
                continue
            pos = f"{chrom}:{shifted_start}-{end}"
            if strand:
                stranded_keys.add(pos + strand)
            else:
                posonly_keys.add(pos)

    return stranded_keys, posonly_keys


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
                "file_path": sample_info.get("file_path"),
                "reference": sample_info.get("reference"),
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

    input_pl = pl.read_csv(broadcast_config.get("file_path"), has_header=False, new_columns=['ID'])

    # Tolerant matching (default behaviour): see _expand_input_keys for details.
    stranded_keys, posonly_keys = _expand_input_keys(input_pl['ID'].to_list())

    reference_col = col(f"{reference}")
    filters = []
    if stranded_keys:
        filters.append(reference_col.is_in(list(stranded_keys)))
    if posonly_keys:
        filters.append(reference_col.str.replace(r'\|[+-]$', '').is_in(list(posonly_keys)))

    if filters:
        filter_expr = filters[0]
        for extra in filters[1:]:
            filter_expr = filter_expr | extra
        hits = lookup_pl.filter(filter_expr)
    else:
        hits = lookup_pl.clear()

    return {sample_name: {db_name: hits}}


def return_collected_results(database_lookup: Collect[Dict[str, Dict[str, pl.DataFrame]]]) -> Dict[str, Dict[str, pl.DataFrame]]:
    """
    Collate sample results into a single dictionary of {sample_name: {db_name: DataFrame}}
    
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


def _reference_hit_sets(hits: pl.DataFrame, reference: str) -> Tuple[set, set]:
    """Build stranded/strand-stripped key sets from a database's already-filtered hits."""
    if hits.is_empty() or reference not in hits.columns:
        return set(), set()
    values = [v for v in hits[reference].to_list() if v is not None]
    stranded = set(values)
    posonly = {re.sub(r'\|[+-]$', '', v) for v in values}
    return stranded, posonly


def _raw_id_matched(raw_id: str, stranded_hit_set: set, posonly_hit_set: set) -> bool:
    """Check whether a raw input ID, with the same +/-1 tolerance as the lookup filter, hit a database."""
    match = _COORD_RE.match(raw_id)
    if match is None:
        return False
    chrom, start, end, strand = match.group(1), int(match.group(2)), match.group(3), match.group(4)
    for shifted_start in (start, start - 1, start + 1):
        if shifted_start < 0:
            continue
        pos = f"{chrom}:{shifted_start}-{end}"
        if strand:
            if pos + strand in stranded_hit_set:
                return True
        elif pos in posonly_hit_set:
            return True
    return False


def _matched_raw_ids(raw_ids: list, hits_by_db: Dict[str, pl.DataFrame], reference: str) -> set:
    """Union raw IDs that hit at least one database's lookup table for a sample."""
    matched: set = set()
    for hits in hits_by_db.values():
        stranded_hit_set, posonly_hit_set = _reference_hit_sets(hits, reference)
        if not stranded_hit_set and not posonly_hit_set:
            continue
        matched.update(
            raw_id for raw_id in raw_ids
            if raw_id not in matched and _raw_id_matched(raw_id, stranded_hit_set, posonly_hit_set)
        )
    return matched


def _write_sample_no_hits(sample_name: str, sample_info: Dict[str, Any], hits_by_db: Dict[str, pl.DataFrame], output_dir: str) -> None:
    """Write unmatched raw IDs for a single sample to no_hits.txt, if any."""
    file_path = sample_info.get("file_path")
    reference = sample_info.get("reference")
    if not file_path or not reference:
        return
    try:
        raw_ids = [ln.strip() for ln in Path(file_path).read_text().splitlines() if ln.strip()]
    except Exception:
        return

    matched = _matched_raw_ids(raw_ids, hits_by_db, reference)
    unmatched = [raw_id for raw_id in raw_ids if raw_id not in matched]
    if not unmatched:
        return

    p = Path(output_dir)
    if p.is_absolute():
        output_path = p / f"{sample_name}_no_hits.txt"
    else:
        output_path = Path(os.getcwd()) / output_dir / f"{sample_name}_no_hits.txt"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(unmatched) + "\n")


def write_no_hits(
    return_collected_results: Dict[str, Dict[str, pl.DataFrame]],
    config: Dict[str, Any],
) -> None:
    """
    Write input circRNAs that matched no database's lookup table to no_hits.txt.

    Depends on the already-collected per-sample hits rather than running its own
    parallel/collect pipeline, since Hamilton's dynamic execution does not support
    requesting two independent Collect() reductions in a single driver.execute() call.
    """
    output_dir = config.get("global_parameters", {}).get("output_dir", "results/")
    for sample_name, sample_info in config.get("samples", {}).items():
        _write_sample_no_hits(sample_name, sample_info, return_collected_results.get(sample_name, {}), output_dir)
