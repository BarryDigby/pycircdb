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

# Initialise the logger
logger = logging.getLogger(__name__)

def stage_mirnas(file_in):

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

        # hardcode for acces i.e miRNA column name 
        id_name = reader.fieldnames[0]

        # Store valid and invalid rows, write if user wants
        rows = []
        mirna_id = []

        # Iterate over dict (rows)
        for row in reader:

            rows.append(row)
            mirna_id.append(row[id_name])

        # strip hsa- from strings
        mirna_id = [x.strip('hsa-') for x in mirna_id]

        # return list of matches 
        matching_mirnas = helpers.check_identifiers("/home/barry/Desktop/pycircdb/test_data/mirna_identifiers.txt", mirna_id)

        # find items in list that are in the set
        in_database = [x for x in mirna_id if x in matching_mirnas]
        not_in_database = [x for x in mirna_id if x not in matching_mirnas]

        if len(in_database) == 0:
            logger.critical("No valid miRNA records found in the input file {}'".format(file_in))
            return {"config": None, "report":None, "sys_exit_code": 1}

        # subset 'rows' list using not_in_database
        # valid rows just means they looked like miRBase IDs, they may not be in the database
        # add dollar to search for end of string
        bad_entries = [d for d in rows if any(re.search(f'{item}$', v) for item in not_in_database for v in d.values())]


        invalid_file = Path(config.output_dir + '/invalid_mirna_records' + ext)
        if len(bad_entries) > 0:
            logger.warning(f"{len(bad_entries)} miRNAs were not found in the database")
            logger.warning(f"Writing invalid miRNA records to {invalid_file}")
            with invalid_file.open("w", newline="") as out_handle:
                writer = csv.DictWriter(out_handle, fieldnames=header, delimiter=delim)
                writer.writeheader()
                for row in bad_entries:
                    writer.writerow(row)

        return matching_mirnas



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
    stage_mirnas()