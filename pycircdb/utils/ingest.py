import os
import csv
import logging
import sys
import re
from collections import defaultdict
from pathlib import Path
import logging
import csv
from . import config, helpers

"""
    Definitions used for ingesting user input to pycircdb
"""

# Initialise the logger
logger = logging.getLogger(__name__)

databases = ['ArrayStar', 'circBank', 'circBase', 'circAtlas']

def check_circrna(file_in):

    # Stage as file and check its extension
    file_in = Path(file_in)
    _, ext = os.path.splitext(file_in)
    assert ext in ['.csv', '.txt', '.tsv'], logger.error("Input circRNA file '{file}' is not a valid text file. Please provide a valid text file with extension '.csv', '.txt', or '.tsv'.")

    # Open file handle - all file reading operations go below. 
    with file_in.open(newline="") as in_handle:
    
        # Stage lines to dict
        dialect = sniff_format(in_handle)
        reader = csv.DictReader(in_handle, dialect=dialect)

        # store for writing corrected files 
        delim = dialect.delimiter
        header = list(reader.fieldnames)

        # hard-code for access
        id_name = reader.fieldnames[0]

        # Stage list to assess what identifiers map to which database
        # Use dict for accessing parquet files later
        # {'key': ['val1', 'val2', 'val3']}.. parquet query accepts list
        database_maps = defaultdict(list)

        # Store valid and invalid rows, write if user wants
        valid_rows = []
        invalid_rows = []

        # Store coordinates seperately. Add to database maps when ref build determined.
        coordinates = []

        # Iterate over dict (rows)
        for row in reader:
            if row[id_name].startswith("ASCRP"):
                valid_rows.append(row)
                database_maps['ArrayStar'].append(row[id_name])
            elif row[id_name].startswith("hsa_circ_") and row[id_name][9].isdigit():
                valid_rows.append(row)
                database_maps["circBase"].append(row[id_name])
            elif row[id_name].startswith("hsa-"):
                valid_rows.append(row)
                database_maps["circAtlas"].append(row[id_name])
            elif re.match(r'hsa_circ[a-zA-Z]\w+', row[id_name]) or row[id_name].startswith('hsa_circ_chr'):
                valid_rows.append(row)
                database_maps["circBank"].append(row[id_name])
            elif re.match(r'^chr([1-9]|1[0-9]|2[0-2]+|X|Y|M):', row[id_name]):
                valid_rows.append(row)
                coordinates.append(row[id_name])
            else:
                invalid_rows.append(row)

        # Let user know which databases their inputs map to
        database_maps = list(set(database_maps))
        if len(database_maps) == 1:
            logger.info("All of the input circRNA identifiers provided map to " + database_maps[0])
        elif len(database_maps) < 1:
            logger.info("User provided a list of circRNAs coordinates as identifiers")
        else:
            logger.info("User provided a mix of circRNA identifiers that map to " + ", ".join(database_maps))

        # Check which reference build the coords match to.
        # clever way to pass arquet file using config here. 
        # User may not know if they are 0-based, 1-based or a mix of both. (CIRI + other tool for e.g)
        if len(coordinates) > 0:
            logger.info("Checking if the input circRNA coordinates belong to Hg19/Hg38.")
            if helpers.query_parquet("test_data/annotations.parquet", "hg19", "in", coordinates):
                logger.info("Input circRNA coordinates belong to Hg19")
            elif helpers.query_parquet("test_data/annotations.parquet", "hg38", "in", coordinates):
                logger.info("Input circRNA coordinates belong to Hg38")
            else:
                logger.info("Input circRNA coordinates do not match. Trying 0-based and 1-based coordinates for your input.")


        # Basic info about inputs
        logger.info("User provided " + str(len(valid_rows)) + " valid circRNA identifiers and " + str(len(invalid_rows)) + " invalid circRNA identifiers.")

        # Write sep files for the user. Should toggle this on/off via param
        valid_out = Path(config.output_dir + "/valid_inputs" + ext)
        invalid_out = Path(config.output_dir + "/invalid_inputs" + ext)

        if len(valid_rows) > 0:
            with valid_out.open(mode="w", newline="") as out_handle:
                writer = csv.DictWriter(out_handle, header, delimiter=delim)
                writer.writeheader()
                for row in valid_rows:
                    writer.writerow(row)

            logger.info( str(len(valid_rows)) + " valid circRNA identifiers have been written to '" + str(valid_out) + "'")

        if len(invalid_rows) > 0:
            with invalid_out.open(mode="w", newline="") as out_handle:
                writer = csv.DictWriter(out_handle, header, delimiter=delim)
                writer.writeheader()
                for row in invalid_rows:
                    writer.writerow(row)

            logger.info( str(len(invalid_rows)) + " invalid circRNA identifiers have been written to '" + str(invalid_out) + "'")





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