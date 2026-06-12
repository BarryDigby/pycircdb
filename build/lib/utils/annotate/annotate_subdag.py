import polars as pl
from polars import col
from hamilton.function_modifiers import config
from typing import Dict, Any, Union, List
from pathlib import Path
import os


PerSampleAnnotation = Dict[str, Union[str, List[str], pl.DataFrame, None]]
AnnotateResult = Dict[str, Dict[str, pl.DataFrame]]


@config.when(db_name='arraystar')
def annotate__arraystar(
    per_sample_annotation: PerSampleAnnotation
) -> AnnotateResult:
    """
    hg19, arraystar
    """
    lookup_pl = per_sample_annotation.get("lookup_hits")
    annotation_pl = pl.read_parquet(per_sample_annotation.get("annotation_parquet_path"))
    annotated_hits = annotation_pl.filter(col("arraystar").is_in(lookup_pl['arraystar']))
    annotated_hits = annotated_hits.join(lookup_pl.select(['arraystar', 'hg38']), on='arraystar', how='left')
    annotated_hits = annotated_hits.select(['arraystar', 'hg19', 'hg38'] + [col for col in annotation_pl.columns if col not in ['arraystar', 'hg19', 'hg38']])
    
    return {per_sample_annotation.get("sample_name"): {'arraystar': annotated_hits}}


@config.when(db_name='circatlas')
def annotate__circatlas(
    per_sample_annotation: PerSampleAnnotation
) -> AnnotateResult:
    """
    circatlas, hg38
    """
    # No database to match against, just return the lookup hits with hg38 coordinates as annotation
    return {per_sample_annotation.get("sample_name"): {'circatlas': per_sample_annotation.get("lookup_hits")}}

@config.when(db_name='circbank')
def annotate__circbank(
    per_sample_annotation: PerSampleAnnotation
) -> AnnotateResult:
    """
    circbank, hg19
    """
    lookup_pl = per_sample_annotation.get("lookup_hits")
    annotation_pl = pl.read_parquet(per_sample_annotation.get("annotation_parquet_path"))
    annotated_hits = annotation_pl.filter(col("circbank").is_in(lookup_pl['circbank']))
    annotated_hits = annotated_hits.join(lookup_pl.select(['circbank', 'hg38']), on='circbank', how='left')
    annotated_hits = annotated_hits.select(['circbank', 'hg19', 'hg38'] + [col for col in annotation_pl.columns if col not in ['circbank', 'hg19', 'hg38']])
    
    return {per_sample_annotation.get("sample_name"): {'circbank': annotated_hits}}


@config.when(db_name='circbase')
def annotate__circbase(
    per_sample_annotation: PerSampleAnnotation
) -> AnnotateResult:
    """
    hg19, name
    """
    lookup_pl = per_sample_annotation.get("lookup_hits")
    annotation_pl = pl.read_parquet(per_sample_annotation.get("annotation_parquet_path"))
    annotated_hits = annotation_pl.filter(col("name").is_in(lookup_pl['circbase']))
    annotated_hits = annotated_hits.join(lookup_pl.select(['circbase', 'hg38']), left_on='name', right_on='circbase', how='left')
    annotated_hits = annotated_hits.rename({'name': 'circbase'})
    annotated_hits = annotated_hits.select(['circbase', 'hg19', 'hg38'] + [column for column in annotation_pl.columns if column not in ['name', 'hg19', 'hg38']])

    return {per_sample_annotation.get("sample_name"): {'circbase': annotated_hits}}


@config.when(db_name='circpedia')
def annotate__circpedia(
    per_sample_annotation: PerSampleAnnotation
) -> AnnotateResult:
    """
    circpedia, hg38
    """
    lookup_pl = per_sample_annotation.get("lookup_hits")
    annotation_pl = pl.read_parquet(per_sample_annotation.get("annotation_parquet_path"))
    annotated_hits = annotation_pl.filter(col("circpedia").is_in(lookup_pl['circpedia']))
    annotated_hits = annotated_hits.join(lookup_pl.select(['circpedia', 'hg19']), on='circpedia', how='left')
    annotated_hits = annotated_hits.select(['circpedia', 'hg19', 'hg38'] + [col for col in annotation_pl.columns if col not in ['circpedia', 'hg19', 'hg38']])

    return {per_sample_annotation.get("sample_name"): {'circpedia': annotated_hits}}


@config.when(db_name='circRNA_DB')
def annotate__circrna_db(
    per_sample_annotation: PerSampleAnnotation
) -> AnnotateResult:
    """
    circRNA_DB, hg19
    """
    lookup_pl = per_sample_annotation.get("lookup_hits")
    annotation_pl = pl.read_parquet(per_sample_annotation.get("annotation_parquet_path"))
    annotated_hits = annotation_pl.filter(col("circRNA_DB").is_in(lookup_pl['circRNA_DB']))
    annotated_hits = annotated_hits.join(lookup_pl.select(['circRNA_DB', 'hg38']), on='circRNA_DB', how='left')
    annotated_hits = annotated_hits.select(['circRNA_DB', 'hg19', 'hg38'] + [col for col in annotation_pl.columns if col not in ['circRNA_DB', 'hg19', 'hg38']])
    
    return {per_sample_annotation.get("sample_name"): {'circRNA_DB': annotated_hits}}


@config.when(db_name='cscd')
def annotate__cscd(
    per_sample_annotation: PerSampleAnnotation
) -> AnnotateResult:
    """
    circRNA (hg38)
    """

    lookup_pl = per_sample_annotation.get("lookup_hits")
    annotation_paths = per_sample_annotation.get("annotation_parquet_path")
    annotated_hits = pl.DataFrame()

    for annotation_path in annotation_paths:
        annotation_pl = pl.read_parquet(annotation_path)
        filtered_annotation = annotation_pl.filter(col("circRNA").is_in(lookup_pl['hg38']))
        if not filtered_annotation.is_empty():
            filtered_annotation = filtered_annotation.join(lookup_pl.select(['hg38', 'hg19']), left_on='circRNA', right_on='hg38', how='left')
            filtered_annotation = filtered_annotation.rename({'circRNA': 'hg38'})
            filtered_annotation = filtered_annotation.select(['hg38', 'hg19'] + [col for col in filtered_annotation.columns if col not in ['hg38', 'hg19']])
            annotated_hits = pl.concat([annotated_hits, filtered_annotation], how='vertical')
        
    return {per_sample_annotation.get("sample_name"): {'cscd': annotated_hits}}


@config.when(db_name='exorbase')
def annotate__exorbase(
    per_sample_annotation: PerSampleAnnotation
) -> AnnotateResult:
    """
    exorbase, hg38
    """
    lookup_pl = per_sample_annotation.get("lookup_hits")
    annotation_pl = pl.read_parquet(per_sample_annotation.get("annotation_parquet_path"))
    annotated_hits = annotation_pl.filter(col("exorbase").is_in(lookup_pl['exorbase']))
    annotated_hits = annotated_hits.join(lookup_pl.select(['exorbase', 'hg19']), on='exorbase', how='left')
    annotated_hits = annotated_hits.select(['exorbase', 'hg19', 'hg38'] + [col for col in annotation_pl.columns if col not in ['exorbase', 'hg19', 'hg38']])

    return {per_sample_annotation.get("sample_name"): {'exorbase': annotated_hits}}


def write_to_output_dir(
    annotate: AnnotateResult,
    per_sample_annotation: PerSampleAnnotation,
) -> None:
    """
    Write annotated hits to output directory as CSV
    Makes sense to do this here to avoid carrying bloat for subdag collection.
    """
    output_dir = per_sample_annotation.get("output_dir")
    sample_name = per_sample_annotation.get("sample_name")
    db_name = per_sample_annotation.get("db_name")

    annotated_hits = annotate.get(sample_name, {}).get(db_name)
    if annotated_hits is None:
        return None
    
    p = Path(output_dir).expanduser()
    if p.is_absolute():
        output_path = p / sample_name
    else:
        cwd_tmp = os.path.join(os.getcwd(), output_dir, sample_name)
        output_path = Path(cwd_tmp)

    output_path.mkdir(parents=True, exist_ok=True)
    annotated_hits.write_csv(output_path / f"{db_name}_hits.txt", separator='\t', include_header=True)

    return None