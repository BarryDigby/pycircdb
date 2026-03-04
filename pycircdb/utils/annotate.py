import os
import logging
from . import config, helpers
import pandas as pd
from collections import defaultdict

# Initialise the logger
logger = logging.getLogger(__name__)

def annotate_circrna(circrna, ann_df):
    """
    Appropriate filtering (strands etc) has been performed in the filter_circrna function.
    Just write the file if filtering and annotation params are provided.
    If they have not been provided, subset the parquet file 'vanilla' and write. 
    """

    if not os.path.exists(f'{config.outdir}/annotate'):
        os.makedirs(f'{config.outdir}/annotate')

    if ann_df is not None:
        ann_df.to_csv(f"{config.outdir}/annotate/circrna_annotations.txt",sep="\t",index=False,na_rep="NA")
    else:

        df = helpers.query_parquet(
            parquet_file="test_data/annotations.parquet",
            key="hg38",
            operator="in",
            value=list(circrna.keys()),
            columns=None,
        ).to_pandas()
    
        df.to_csv(f"{config.outdir}/annotate/circrna_annotations.txt",sep="\t",index=False,na_rep="NA")

    logger.info(f"Annotation file created: {config.outdir}/annotate/circrna_annotations.txt")


if __name__ == "__main__":
    annotate_circrna()
