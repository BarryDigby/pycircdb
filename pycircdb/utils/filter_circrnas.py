import os
import logging
from . import config, helpers
import pandas as pd
from collections import defaultdict

"""
    Filter the users input list of valid circRNA identifiers.
    Takes as input the circrna dictionary of {'hg38': 'original user input'}
    parse_identifiers is used to separate strands from positional IDs (database does not contain strands in ID)
    filter_by_strand creates tmp columns ID + strand which are used to filter the annotation file
"""

# Initialise the logger
logger = logging.getLogger(__name__)

def filter_circrna(circrna_dict, algorithms=None, set_logic=None):

    # sanity checks
    print("printing algorithms")
    print(algorithms)
    print("printing set_logic")
    print(set_logic)

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
        dataframe_collection.append(filter_by_strand(identifier_strand))

    master = (
        pd.concat(dataframe_collection)
        .sort_values(["hg19", "hg38"])
        .reset_index(drop=True)
    )

    if set_logic == "OR":
        subset_df = master[master["algorithm"].str.contains("|".join(algorithms), case=False, na=False)]
    else:
        mask = pd.Series([all(word in item for word in algorithms)for item in master["algorithm"]])
        subset_df = master[mask]

    circrna_dict = {k: v for k, v in circrna_dict.items() if k in subset_df["hg38"].tolist()}

    return circrna_dict, subset_df


def filter_by_strand(identifier_strand):
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
    filter_circrna()
