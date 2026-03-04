#!/usr/bin/env python

import logging
from . import config, helpers
import pandas as pd
import dask.dataframe as dd
from distributed import Client
from pathlib import Path

"""
Inputs:
    circrna:
        {'hg38 identifier': 'original user input'}
    
    mirna/mrna/rbp:
        ['miRNA', 'miRNA',...]
"""

# Initialise the logger
logger = logging.getLogger(__name__)

def initialise_network(circrna=None, mirna=None, gene=None, rbp=None, workers=None):

    # We need to assess what inputs the user has provided
    # Direction of network building dictated by this.
    # if circrna is not None, set boolean flag to True
    # if mirna is not None, set boolean flag to True
    circrna_flag = True if circrna is not None else False
    mirna_flag = True if mirna is not None else False
    gene_flag = True if gene is not None else False
    rbp_flag = True if rbp is not None else False
    
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
    initialise_network()