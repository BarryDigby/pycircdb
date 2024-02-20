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


    # build simple network - circRNA as input
    # testing dask here...
    circrna_mirna(circrna, mirna)
    


def circrna_mirna(circrna, mirna):

    client = Client(n_workers=4, threads_per_worker=1)

    client.close()




if __name__ == "__main__":
    tester()