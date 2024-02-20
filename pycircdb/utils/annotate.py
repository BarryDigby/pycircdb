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


def annotate_circRNAs(annotate, circrna_dict, algorithms=None, set_logic=None):
    # {'hg38_id':'original_input_id'...}
    # take care to note strand information to subset results (prior to CLI param subsets)
    # create coordinate DF, subset by strand if required, add to other annotation dbs. 

    contains_strand = defaultdict(list)
    no_strand = []
    collect_df = []

    for hg38, original_id in circrna_dict.items():
        # check if strands are present if positonal identifiers used
        if ":+" in original_id or ":-" in original_id:
            contains_strand[hg38].append(original_id.split(":")[2])
        else:
            no_strand.append(hg38)

    # this will not need to be run if DB ID's are input
    if len(contains_strand) > 0:
        strand_df = helpers.query_parquet(
            parquet_file="test_data/annotations.parquet",
            key='hg38',
            operator="in",
            value=list(contains_strand.keys()),
            columns=None,
            ).to_pandas()

        # create hg19/hg38 + strand column
        strand_df["tmp"] = strand_df["hg19"] + ":" + strand_df["strand"]
        strand_df["tmp2"] = strand_df["hg38"] + ":" + strand_df["strand"]

        # 'expand' dictionary key vals to string list to check against tmp columns
        # only append ":{strand}" to position if strand is present
        lst = []
        for position, strand in contains_strand.items():
            if strand[0] in ['+', '-']:
                lst.append(f"{position}:{strand[0]}")
            else:
                lst.append(f"{position}")

        strand_df = strand_df.loc[strand_df["tmp"].isin(lst) | strand_df["tmp2"].isin(lst)]

        # drop the temporary columns
        strand_df.drop(columns=["tmp", "tmp2"], inplace=True)

        collect_df.append(strand_df)

    # now grab the other annotations and concat them 
    df =  helpers.query_parquet(
        parquet_file="test_data/annotations.parquet",
        key='hg38',
        operator="in",
        value=no_strand,
        columns=None,
        ).to_pandas()
    
    collect_df = [df] + collect_df

    master = pd.concat(collect_df).sort_values(['hg19', 'hg38']).reset_index(drop=True)

    if annotate:
        logger.info(f"Annotation file created: {config.output_dir}/circrna_annotations.txt")
        master.to_csv(f"{config.output_dir}/circrna_annotations.txt",sep="\t",index=False,na_rep="NA")

    # Now create subset and export this to downstream modules. Write it for the user too.
    if algorithms is not None and set_logic is None:
        mask = pd.Series([all(word in item for word in algorithms) for item in master['algorithm']])
        subset_df = master[mask]
        if annotate:
            subset_df.to_csv(f"{config.output_dir}/filtered_circrna_annotations_subset.txt", sep="\t", index=False, na_rep="NA")
    elif algorithms is not None and set_logic is not None:
        if set_logic == "OR":
            subset_df = master[master['algorithm'].str.contains('|'.join(algorithms), case=False, na=False)]
            if annotate:
                subset_df.to_csv(f"{config.output_dir}/filtered_circrna_annotations_subset.txt", sep="\t", index=False, na_rep="NA")
        elif set_logic == "AND":
            mask = pd.Series([all(word in item for word in algorithms) for item in master['algorithm']])
            subset_df = master[mask]
            if annotate:
                subset_df.to_csv(f"{config.output_dir}/filtered_circrna_annotations_subset.txt", sep="\t", index=False, na_rep="NA")

    # subset the original input if algorithms were appplied.
    if algorithms is not None:
        circrna_dict = {k: v for k, v in circrna_dict.items() if k in subset_df['hg38'].tolist()}

    return circrna_dict


if __name__ == "__main__":
    annotate_circRNAs()
