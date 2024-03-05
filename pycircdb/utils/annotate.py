import os
import logging
from . import config, helpers
import pandas as pd
from collections import defaultdict

"""
    Takes as input a dictionary.
    Subset the annotation parquet file to return annotation file for user
    Schema:
    ['hg19', 'hg38', 'circBase_ID', 'circAtlas_ID', 'circBank_ID', 'Arraystar_ID', 'HGNC', 'strand', 'spliced_length', 'circBase_type', 'CSCD_type', 'algorithm']
"""

# Initialise the logger
logger = logging.getLogger(__name__)


def annotate_filter_circrnas(annotate, circrna_dict, algorithms=None, set_logic=None):

    identifier_strand = defaultdict(list)
    identifier = []
    dataframe_collection = []

    identifier_strand, identifier = parse_identifiers(
        circrna_dict, identifier_strand, identifier
    )

    df = helpers.query_parquet(
        parquet_file="test_data/annotations.parquet",
        key="hg38",
        operator="in",
        value=identifier,
        columns=None,
    ).to_pandas()

    dataframe_collection.append(df)

    if len(identifier_strand) > 0:
        dataframe_collection.append(annotate_filter_by_strand(identifier_strand))

    master = (
        pd.concat(dataframe_collection)
        .sort_values(["hg19", "hg38"])
        .reset_index(drop=True)
    )

    if algorithms is not None:

        if set_logic == "OR":
            subset_df = master[master["algorithm"].str.contains("|".join(algorithms), case=False, na=False)]
        else:
            mask = pd.Series([all(word in item for word in algorithms)for item in master["algorithm"]])
            subset_df = master[mask]

        circrna_dict = {k: v for k, v in circrna_dict.items() if k in subset_df["hg38"].tolist()}

    if annotate:
        if not os.path.exists(f'{config.outdir}/annotate'):
            os.makedirs(f'{config.outdir}/annotate')
        
        logger.info(f"Annotation file created: {config.outdir}/annotate/circrna_annotations.txt")
        
        master.to_csv(f"{config.outdir}/annotate/circrna_annotations.txt",sep="\t",index=False,na_rep="NA")

        if algorithms is not None:
            logger.info(f"Filtered annotation file created: {config.outdir}/annotate/circrna_annotations.txt")
            subset_df.to_csv(f"{config.outdir}/annotate/filtered_circrna_annotations.txt",sep="\t",index=False,na_rep="NA")

    return circrna_dict


def annotate_filter_by_strand(identifier_strand):
    identifiers_with_strand_df = helpers.query_parquet(
        parquet_file="test_data/annotations.parquet",
        key="hg38",
        operator="in",
        value=list(identifier_strand.keys()),
        columns=None,
    ).to_pandas()

    # create hg19/hg38 + strand column
    identifiers_with_strand_df["tmp"] = (
        identifiers_with_strand_df["hg19"] + ":" + identifiers_with_strand_df["strand"]
    )
    identifiers_with_strand_df["tmp2"] = (
        identifiers_with_strand_df["hg38"] + ":" + identifiers_with_strand_df["strand"]
    )

    # 'expand' dictionary key vals to string list to check against tmp columns
    # only append ":{strand}" to position if strand is present
    lst = []
    for position, strand in identifier_strand.items():
        if strand[0] in ["+", "-"]:
            lst.append(f"{position}:{strand[0]}")
        else:
            lst.append(f"{position}")

    identifiers_with_strand_df = identifiers_with_strand_df.loc[
        identifiers_with_strand_df["tmp"].isin(lst)
        | identifiers_with_strand_df["tmp2"].isin(lst)
    ]

    # drop the temporary columns
    identifiers_with_strand_df.drop(columns=["tmp", "tmp2"], inplace=True)

    return identifiers_with_strand_df


def parse_identifiers(
    circrna_dict, positional_strand_information, non_positional_identifiers
):
    for hg38, original_id in circrna_dict.items():
        if ":+" in original_id or ":-" in original_id:
            positional_strand_information[hg38].append(original_id.split(":")[2])
        else:
            non_positional_identifiers.append(hg38)

    return positional_strand_information, non_positional_identifiers


if __name__ == "__main__":
    annotate_filter_circrnas()
