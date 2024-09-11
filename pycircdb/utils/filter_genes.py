import os
import logging
from . import config, helpers
import dask.dataframe as dd
from distributed import Client 
import pandas as pd
from collections import defaultdict

"""
    Filter miRNAs for user, return a list.
    Need to use Dask (read parquet, not pyarrow) here to dynamically apply filtering functions on the circrna-mirna parquet files.
    miRNA, HGNC, database, experiment, support - strings
"""

# Initialise the logger
logger = logging.getLogger(__name__)

def filter_gene(gene_list,
                databases=None,
                set_logic=None,
                workers=None):

    # sanity checks
    #print("printing gene databases") 
    #print(databases)
    #print("printing set_logic")
    #print(set_logic)
#
    client = Client(n_workers=workers, threads_per_worker=1)

    gene_filtering_variables = [databases]
    if any(var is not None for var in gene_filtering_variables):

        # Do not ingest the entire parquet file! 
        columns_dict = {
            databases: 'database',
        }

        filter_function_dict = {
            'HGNC': (filter_by_gene, (gene_list,)),
            'database': (filter_gene_by_database, (databases, set_logic))
        }

        cols = [ 'miRNA','HGNC'] + [col for var, col in columns_dict.items() if var is not None]

        # use try, these could fail and client.close() would not be called
        try:
            ddf = dd.read_parquet("test_data/mirna_mrna.parquet", engine="pyarrow")
            for col, (func, args) in filter_function_dict.items():
                if col in cols:
                    ddf = ddf.map_partitions(func, *args, meta=ddf)
            ddf = ddf.compute()
            ddf.sort_values(['miRNA', 'HGNC']).reset_index(drop=True)
            ddf.to_csv("results/filtered_genes.txt", sep="\t", index=False)
            client.close()
            return ddf
        except:
            client.close()
            return None

    else:
        # we do nothing, just assign original input gene_list to filtered_genes.
        client.close()
        return gene_list


def filter_by_gene(df, gene_list):
    """
    Filter the dataframe by the HGNC IDs provided by the user.
    Only run if other filtering DB param provided. 
    """
    subset_df = df[df['HGNC'].isin(gene_list)]
    return subset_df

def filter_gene_by_database(df, databases, set_logic):
    """
    Filter the dataframe by the databases provided by the user.
    """
    if set_logic == "OR":
        subset_df = df[df["database"].str.contains("|".join(databases), case=False, na=False)]
    else:
        df = df.reset_index(drop=True)
        mask = pd.Series([all(word in item for word in databases) for item in df["database"]])
        subset_df = df[mask]
    return subset_df

if __name__ == "__main__":
    filter_gene()
