import os
import csv
import logging
import sys
import re
from collections import defaultdict
from pathlib import Path
import logging
import csv
from .. import config, helpers
import pyarrow.parquet as pq
import pandas as pd
import dask.dataframe as dd
import dask.distributed
import multiprocessing as mp

# Initialise the logger
logger = logging.getLogger(__name__)

def process_chunk(chunk, list):
    return [d for d in chunk if any(item in v for item in list for v in d.values())]

def stage_mrnas(file_in):

    # Stage as file and check its extension
    file_in = Path(file_in)
    _, ext = os.path.splitext(file_in)
    assert ext in ['.csv', '.txt', '.tsv'], logger.error("Input miRNA file '{file}' is not a valid text file. Please provide a valid text file with extension '.csv', '.txt', or '.tsv'.")

    # Open file handle - all file reading operations go below. 
    with file_in.open(newline="") as in_handle:
    
        # Stage lines to dict
        dialect = sniff_format(in_handle)
        reader = csv.DictReader(in_handle, dialect=dialect)

        # store for writing corrected files 
        delim = dialect.delimiter
        header = list(reader.fieldnames)

        # hardcode for acces i.e mRNA column name 
        id_name = reader.fieldnames[0]

        # Store rows, identifiers
        rows = []
        mrna_id = []


        # Iterate over dict (rows)
        for row in reader:

            rows.append(row)
            mrna_id.append(row[id_name])
        
        # return list of matches
        # test data contains a lot of miRNA, pseudogenes and lncRNAs! A good example of non-database matches for circRNAs
        matching_mrnas = helpers.check_identifiers("/home/barry/Desktop/pycircdb/test_data/mrna_identifiers.txt", mrna_id)

        # Convert matching_mrnas to a set for faster lookups
        matching_mrnas_set = set(matching_mrnas)

        # Use a single loop to populate in_database and not_in_database
        in_database = []
        not_in_database = []
        for x in mrna_id:
            if x in matching_mrnas_set:
                in_database.append(x)
            else:
                not_in_database.append(x)

        if len(in_database) == 0:
            logger.critical("No valid mRNA records found in the input file {}'".format(file_in))
            return {"config": None, "report":None, "sys_exit_code": 1}

        # Split the data into chunks
        chunk_size = len(rows) // mp.cpu_count()
        chunks = [rows[i:i + chunk_size] for i in range(0, len(rows), chunk_size)]

        # Create a pool of workers
        with mp.Pool() as pool:
            results = pool.starmap(process_chunk, [(chunk, not_in_database) for chunk in chunks])

        bad_entries = [item for sublist in results for item in sublist]

        #bad_entries = [d for d in rows if any(item in v for item in not_in_database for v in d.values())] 
        #bad_entries = [d for d in rows if any(re.search(f'{item}$', v) for item in not_in_database for v in d.values())]

        invalid_file = Path(config.output_dir + '/invalid_mrna_records' + ext)
        if len(bad_entries) > 0:
            logger.warning(f"{len(bad_entries)} mRNAs were not found in the database")
            logger.warning(f"Writing invalid mRNA records to {invalid_file}")
            with invalid_file.open("w", newline="") as out_handle:
                writer = csv.DictWriter(out_handle, fieldnames=header, delimiter=delim)
                writer.writeheader()
                for row in bad_entries:
                    writer.writerow(row)

        return matching_mrnas



def read_head(handle, num_lines=10):
    """
    Read n lines from current position in file
    """
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)

def sniff_format(handle):
    """
    Detect the tabular format.
    Does not work with single column file
    work around: https://bugs.python.org/issue2078:
    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    if not sniffer.has_header(peek):
        logger.critical(f"The input circRNA file does not have a header")
        sys.exit(1)
    dialect = sniffer.sniff(peek)
    if dialect.delimiter not in [",", "\t"]:
        logger.debug(f"Single column file provided by user. Default to tab delimiter")
        dialect.delimiter = "\t"

    return dialect


if __name__ == "__main__":
    stage_mrnas()