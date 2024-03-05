#!/usr/bin/env python

import logging
from . import config, helpers
import pandas as pd
import dask.dataframe as dd
from distributed import Client
from pathlib import Path

"""
    go back to circrna ingest and make a key value pair hg38 > original input?
    or do you just go forward with hg38 , any need to preserve the original probe ID? 
    yes there is actually if they want to keep ASC etc.. 
"""

# Initialise the logger
logger = logging.getLogger(__name__)


def tester(circrna=None, mirna=None):


    circrna_mirna(circrna)


def filter_matching_circrnas(dask_df, circrna_ids):
    """
    Filter the circrna ids from the dask dataframe
    """
    dask_df = dask_df[dask_df['hg38'].isin(circrna_ids)]
    return dask_df



def circrna_mirna(circrna):

    client = Client(n_workers=4, threads_per_worker=1)

    try:
        ddf = dd.read_parquet("test_data/mirna_chrs/*.parquet", engine="pyarrow", columns=['hg38', 'miRNA'])

        ddf = ddf.map_partitions(filter_matching_circrnas, list(circrna.keys())).compute()
        print(len(ddf))

        ddf.to_csv("results/filtered_mirna.csv", index=False)


        client.close()
    
    except:
        client.close()




if __name__ == "__main__":
    tester()