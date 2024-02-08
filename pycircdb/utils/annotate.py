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

def annotate_circRNAs(circrna_dict):

    collect_df = []

    # Iterate over the key-value pairs in the dictionary
    for database, identifier in circrna_dict.items():

        # check if strands are present if positonal identifiers used
        if ':+' in identifier or ':-' in identifier:
            position_strand = defaultdict(list)
            for item in identifier:
                if len(item.split(':')) == 2:
                    position = f'{item.split(":")[0]}:{item.split(":")[1]}'
                    position_strand[position].append('')
                else:
                    position = f'{item.split(":")[0]}:{item.split(":")[1]}'
                    position_strand[position].append(item.split(":")[2])

            # note user may have supplied coordinate with two strands (circRNA_finder)
            positions = list(position_strand.keys())

            #subset the parquet file
            df = helpers.query_parquet(
            parquet_file='test_data/annotations.parquet',
            key=database,
            operator='in',
            value=positions,
            columns=None).to_pandas()

            # create hg19/hg38 + strand column
            df['tmp'] = df['hg19'] + ':' + df['strand']
            df['tmp2'] = df['hg38'] + ':' + df['strand']

            # 'expand' dictionary key vals to string list to check against tmp columns
            lst = [f'{position}:{strand}' for position, values in position_strand.items() for strand in values]
            df = df.loc[df['tmp'].isin(lst) | df['tmp2'].isin(lst)]

            # drop the temporary columns
            df.drop(columns=['tmp', 'tmp2'], inplace=True)

            collect_df.append(df)


        #subset the parquet file
        df = helpers.query_parquet(
            parquet_file='test_data/annotations.parquet',
            key=database,
            operator='in',
            value=identifier,
            columns=None).to_pandas()

        collect_df.append(df)

    # Concatenate the dataframes, output for user.
    # Makes sense to drop duplicates entries here 
    master = pd.concat(collect_df).sort_values(['hg19', 'hg38']).drop_duplicates().reset_index(drop=True)

    logger.info(f'Annotation file created: {config.output_dir}/circrna_annotations.txt')

    master.to_csv(f'{config.output_dir}/circrna_annotations.txt', sep="\t", index=False, na_rep="NA")

if __name__ == "__main__":
    annotate_circRNAs()