#!/usr/bin/env python


import logging
from . import config, helpers
import pandas as pd

"""
    Set of functions for constructing a circRNA-miRNA-mRNA network

    1. Basic - just circRNA inputs. Network can be refined by circRNA algorithm confidence, 
        miRNA prediction tool confidence, and miRNA-mRNA interaction confidence.

    2. User provides input miRNA or both mRNA files which will be used to subset the network. 


    This will get complicated when user has different inputs i.e starting at mRNA..
    Just start coding it up and logic will become more apparent. 
"""

# Initialise the logger
logger = logging.getLogger(__name__)


def stage_inputs(circrna=None,
                mirna=None, 
                mrna=None, 
                circ_algorithm=None, 
                mirna_algorithm=None, 
                mrna_algorithm=None):


    # User has provided circRNA file and wants to build network without filters
    if circrna and all(var is None for var in [mirna, mrna, circ_algorithm, mirna_algorithm, mrna_algorithm]):
        print(circrna)



if __name__ == "__main__":
    stage_inputs()