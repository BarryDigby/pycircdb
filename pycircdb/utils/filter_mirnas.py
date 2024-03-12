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
    hg38/miRNA/algorithm/MRE_start/MRE_end/site_type - string, score/MFE - double
    Decide what the output will be
"""

# Initialise the logger
logger = logging.getLogger(__name__)

def filter_mirna(mirna_list,
                algorithms=None,
                set_logic=None,
                type=None,
                mfe=None,
                score=None,
                workers=None):

    # sanity checks
    print("printing algorithms") 
    print(algorithms)
    print("printing set_logic")
    print(set_logic)
    print("printing type")
    print(type)
    print("printing mfe")
    print(mfe)
    print("printing score")
    print(score)

    mirna_filtering_variables = [algorithms, type, mfe, score]
    if any(var is not None for var in mirna_filtering_variables):

        client = Client(n_workers=workers, threads_per_worker=1)

        # Do not ingest the entire parquet file! 
        columns_dict = {
            algorithms: 'algorithm',
            type: 'site_type',
            mfe: 'MFE',
            score: 'score'
        }

        filter_function_dict = {
            'miRNA': (filter_by_mirnas, (mirna_list,)),
            'algorithm': (filter_mirna_by_algorithm, (algorithms, set_logic)),
            'site_type': (filter_mirna_by_type, (type,)),
            'MFE': (filter_mirna_by_mfe, (mfe,)),
            'score': (filter_mirna_by_score, (score,))
        }

        cols = [ 'hg38','miRNA'] + [col for var, col in columns_dict.items() if var is not None]

        # use try, these could fail and client.close() would not be called
        try:
            # Read all columns from the parquet files
            ddf = dd.read_parquet("test_data/mirna_chrs/*.parquet", engine="pyarrow")
            for col, (func, args) in filter_function_dict.items():
                if col in cols:
                    # Apply the filters only on the necessary columns
                    ddf = ddf.map_partitions(func, *args, meta=ddf)
            ddf = ddf.compute()
            ddf.sort_values(['hg38', 'miRNA']).reset_index(drop=True)
            ddf.to_csv("results/filtered_mirna.txt", sep="\t", index=False)
            client.close()
            return ddf
        except:
            client.close()
            return None
    
    else:
        # we do nothing, just assign original input mirna_list to filtered_mirnas.
        return mirna_list


def filter_by_mirnas(df, mirna_list):
    """
    Filter the dataframe by the miRNAs provided by the user.
    This only needs to be run if other filtering params are provided to reduce comp overhead.
    """
    subset_df = df[df['miRNA'].isin(mirna_list)]
    return subset_df

def filter_mirna_by_algorithm(df, algorithms, set_logic):
    """
    Filter the dataframe by the algorithm provided by the user.
    NaNs exist in algorithm - revise and remove these from DB files as they contain no useful information. 
    """
    if set_logic == "OR":
        subset_df = df[df["algorithm"].str.contains("|".join(algorithms), case=False, na=False)]
    else:
        df = df.reset_index(drop=True)
        #mask = pd.Series([all(word in item for word in algorithms) if pd.notna(item) else False for item in df["algorithm"]])
        # test without nadrop 
        mask = pd.Series([all(word in item for word in algorithms) for item in df["algorithm"]])
        subset_df = df[mask]
    return subset_df

def filter_mirna_by_type(df, type):
    """
    Filter the dataframe by the site type provided by the user.
    """
    type = list(type)
    subset_df = df[df['site_type'].isin(type)]
    return subset_df

def filter_mirna_by_mfe(df, mfe):
    """
    Filter the dataframe by the MFE provided by the user.
    MFE is negative, -62.0 - -0.41 (kcal/mol)
    Use less than to filter
    """
    subset_df = df[df['MFE'] <= mfe]
    return subset_df

def filter_mirna_by_score(df, score):
    """
    Filter the dataframe by the score provided by the user.
    Score is positive, 140.0 - 220.0
    Use greater than to filter
    """
    subset_df = df[df['score'] >= score]
    return subset_df

if __name__ == "__main__":
    filter_mirna()
