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
import pyarrow.parquet as pq
import pandas as pd

"""
    Definitions used for ingesting user input to pycircdb
"""

# Initialise the logger
logger = logging.getLogger(__name__)

databases = ['ArrayStar', 'circBank', 'circBase', 'circAtlas', 'no database (coordinates)']

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
                database_maps["no database (coordinates)"].append(row[id_name])
            else:
                invalid_rows.append(row)

        # Let user know which databases their inputs map to
        if len(list(set(database_maps))) == 1:
            logger.info("All of the input circRNA identifiers provided map to " + list(set(database_maps))[0])
        else:
            logger.info("User provided a mix of circRNA identifiers that map to " + ", ".join(list(set(database_maps))))

        # Check which reference build the coords match to.
        # User may not know if they are 0-based, 1-based or a mix of both. (CIRI + other tool for e.g)
        replacements = defaultdict(list)
        if len(coordinates) > 0:
            logger.info("Checking if the input circRNA coordinates belong to Hg19")

            hg19_matches = helpers.query_parquet("test_data/annotations.parquet", "hg19", "in", coordinates, ['hg19']).to_pandas()
            hg19_matches = hg19_matches['hg19'].to_list()
            hg19_matches = sorted(list(set(hg19_matches)))

            logger.info("Checking if the input circRNA coordinates belong to Hg38")
            hg38_matches = helpers.query_parquet("test_data/annotations.parquet", "hg38", "in", coordinates, ['hg38']).to_pandas()
            hg38_matches = hg38_matches['hg38'].to_list()
            hg38_matches = sorted(list(set(hg38_matches)))

            coordinates = sorted(list(set(coordinates)))

            # Best case scenario, inputs are sound:
            if all(item in hg19_matches + hg38_matches for item in coordinates):
                if hg19_matches is not None:
                    logger.info("Input circRNA coordinates belong to Hg19")
                    database_maps['hg19'] = hg19_matches
                    del database_maps['no database (coordinates)']
                elif hg38_matches is not None:
                    logger.info("Input circRNA coordinates belong to Hg38")
                    database_maps['hg38'] = hg38_matches
                    del database_maps['no database (coordinates)']
            
            # Fix coordinates and attempt to map for user:
            elif not all(item in hg19_matches + hg38_matches for item in coordinates):

                # Stage the database map here, then extract problem coordinates
                database_maps['hg19'] = hg19_matches
                database_maps['hg38'] = hg38_matches
                del database_maps['no database (coordinates)']

                # Extract problem coordinates
                problem_coordinates = [item for item in coordinates if item not in hg19_matches + hg38_matches]
                logger.info(str(len(problem_coordinates)) + " input circRNA coordinates do not match Hg19 or Hg38. Trying 0-based and 1-based coordinates for your input.")
                
                hg19_counter = 0
                hg38_counter = 0
                for coordinate in problem_coordinates:

                    zero_based = adjust_coordinates([coordinate], 'subtract')
                    one_based = adjust_coordinates([coordinate], 'add')

                    zero_19 = helpers.query_parquet("test_data/annotations.parquet", "hg19", "in", zero_based, ['hg19']).to_pandas()['hg19'].drop_duplicates().tolist()
                    one_19 = helpers.query_parquet("test_data/annotations.parquet", "hg19", "in", one_based, ['hg19']).to_pandas()['hg19'].drop_duplicates().tolist()

                    if len(zero_19 + one_19) > 0:
                        hg19_counter += 1
                
                    zero_38 = helpers.query_parquet("test_data/annotations.parquet", "hg38", "in", zero_based, ['hg38']).to_pandas()['hg38'].drop_duplicates().tolist()
                    one_38 = helpers.query_parquet("test_data/annotations.parquet", "hg38", "in", one_based, ['hg38']).to_pandas()['hg38'].drop_duplicates().tolist()

                    if len(zero_38 + one_38) > 0:
                        hg38_counter += 1

                    # Corrected cooridinate in dict, append to database maps. 
                    var_dict = {'hg19': zero_19 + one_19, 'hg38': zero_38 + one_38}
                    if len(var_dict['hg19']) > 0:
                        database_maps['hg19'].append(var_dict['hg19'][0])
                    else:
                        database_maps['hg38'].append(var_dict['hg38'][0])

                    # Isolate corrected coordinate, use to re-write user input file
                    # {'old_coord': ['new_coord']}
                    corrected_coordinate = [list for list in var_dict.values() if list]
                    replacements[coordinate] = corrected_coordinate[0][0]


                if hg19_counter > 0 and hg38_counter == 0:
                    logger.info("circRNA coordinates corrected and mapped to Hg19")
                elif hg19_counter == 0 and hg38_counter > 0:
                    logger.info("circRNA coordinates corrected and mapped to Hg38")
                else:
                    logger.info("circRNA coordinates corrected and mapped to both Hg19 and Hg38")

                # Now perform "surgery" on the input file to inject the correct coordinates for the user
                # must remain within this scope (arraystar probes wont use replacement variable)
                for row in valid_rows:
                    # key = orig, value = corrected
                    for key, values in replacements.items():
                        if row[id_name] == key:
                            row[id_name] = values


        # Basic info about inputs
        logger.info("User provided " + str(len(valid_rows)) + " valid circRNA identifiers and " + str(len(invalid_rows)) + " invalid circRNA identifiers.")

        # Write sep files for the user. Should toggle this on/off via param
        valid_out = Path(config.output_dir + "/valid_inputs" + ext)
        invalid_out = Path(config.output_dir + "/invalid_inputs" + ext)

        if len(valid_rows) > 0 and len(replacements) > 0:
            with valid_out.open(mode="w", newline="") as out_handle:
                writer = csv.DictWriter(out_handle, header, delimiter=delim)
                writer.writeheader()
                for row in valid_rows:
                    writer.writerow(row)

            logger.info("Original input (" + str(len(coordinates) - len(replacements)) +") and corrected circRNA identifiers (" + str(len(replacements)) + ") have been written to '" + str(valid_out) + "'")

        if len(invalid_rows) > 0:
            with invalid_out.open(mode="w", newline="") as out_handle:
                writer = csv.DictWriter(out_handle, header, delimiter=delim)
                writer.writeheader()
                for row in invalid_rows:
                    writer.writerow(row)

            logger.info( str(len(invalid_rows)) + " invalid circRNA identifiers have been written to '" + str(invalid_out) + "'")

        
        # Return the database map. Perform 'surgery' on it upstream within the coordinate scope
        # Remove a key if its empty
        database_maps = dict((k, v) for k, v in database_maps.items() if v)
        
        return database_maps



def adjust_coordinates(coord_list, operation):
    adjusted_coords = []
    for coord in coord_list:
        chr_part, num_part = coord.split(':')
        start, end = num_part.split('-')
        if operation == 'add':
            adjusted_start = str(int(start) + 1)
        elif operation == 'subtract':
            adjusted_start = str(int(start) - 1)
        else:
            raise ValueError("Invalid operation. Choose 'add' or 'subtract'.")
        adjusted_coord = ':'.join([chr_part, '-'.join([adjusted_start, end])])
        adjusted_coords.append(adjusted_coord)
    return adjusted_coords


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