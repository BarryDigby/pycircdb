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
import pandas as pd
import multiprocessing as mp

# Initialise the logger
logger = logging.getLogger(__name__)

def stage_inputs(circrna_file=None, mirna_file=None, gene_file=None, rbp_file=None, workers=None):
    """
    Script entry point, utilises functions below
    to remove entries from user that != databases.
    Makes effort to correct when circRNAs provided.

    Outputs:

    circna: dict
        Dictionary containing {'hg38_id':'original_input_id'...}
        Used to build networks and write data using original user input.
    mirna/mrna/rbp: list
        A list containing identifiers present in the databases and thus
        valid for downstream network construction.
    """

    circrna_dict = mirna_list = gene_list = rbp_list = None

    if circrna_file is not None:
        logger.info("Ingesting circRNA file...")
        circrna_dict = ingest_circrna(circrna_file)

    if mirna_file is not None:
        logger.info('Ingesting miRNA file...')
        mirna_list = ingest_others(mirna_file, workers, "miRNA")

    if gene_file is not None:
        logger.info('Ingesting genes...')
        gene_list = ingest_others(gene_file, workers, "gene")

    if rbp_file is not None:
        logger.info('Ingesting RBPs...')
        rbp_list = ingest_others(rbp_file, workers, "RBP")

    return circrna_dict, mirna_list, gene_list, rbp_list

def ingest_circrna(file_in):
    """
    Frankenstein code, forgive me... 
    """

    file_path = Path(file_in)
    _, ext = os.path.splitext(file_path)

    with file_path.open(newline="") as in_handle:

        file_dialect = sniff_format(in_handle)
        file_reader = csv.DictReader(in_handle, dialect=file_dialect)
        delim = file_dialect.delimiter
        header = list(file_reader.fieldnames)

        identifier_colname = file_reader.fieldnames[0]

        database_maps = defaultdict(list)
        valid_rows = []
        invalid_rows = []
        coordinates = []

        for row in file_reader:

            if row[identifier_colname].startswith("ASCRP"):
                valid_rows.append(row)
                database_maps["ArrayStar"].append(row[identifier_colname])
            elif row[identifier_colname].startswith("hsa_circ_") and row[identifier_colname][9].isdigit():
                valid_rows.append(row)
                database_maps["circBase"].append(row[identifier_colname])
            elif row[identifier_colname].startswith("hsa-"):
                valid_rows.append(row)
                database_maps["circAtlas"].append(row[identifier_colname])
            elif re.match(r"hsa_circ[a-zA-Z]\w+", row[identifier_colname]) or row[identifier_colname].startswith("hsa_circ_chr"):
                valid_rows.append(row)
                database_maps["circBank"].append(row[identifier_colname])
            elif re.match(r"^chr([1-9]|1[0-9]|2[0-2]+|X|Y|M):", row[identifier_colname]):
                valid_rows.append(row)
                coordinates.append(row[identifier_colname])
                database_maps["no database (coordinates)"].append(row[identifier_colname])
            else:
                invalid_rows.append(row)

        # Check which reference build the coords match to.
        # User may not know if they are 0-based, 1-based or a mix of both. (CIRI + other tool for e.g)
        replacements = defaultdict(list)
        if len(coordinates) > 0:
            logger.debug("Checking if the input circRNA coordinates belong to Hg19")

            # Remove the strand info when performing DB matches.
            # but store in a dict as we need it later.
            strand_info = defaultdict(list)
            for i, coordinate in enumerate(coordinates):
                coordinate_parts = coordinate.split(":")
                if len(coordinate_parts) == 3:
                    chr_part, num_part, strand = coordinate_parts
                    coordinate_less_strand = ":".join([chr_part, num_part])
                    coordinates[i] = coordinate_less_strand
                    strand_info[coordinate_less_strand] = strand

            # Use a set() so duplicate positions on + & - strand are not a problem.
            hg19_matches = helpers.query_parquet(
                "test_data/annotations.parquet", "hg19", "in", coordinates, ["hg19"]
            ).to_pandas()
            hg19_matches = hg19_matches["hg19"].to_list()
            hg19_matches = sorted(list(set(hg19_matches)))

            logger.debug("Checking if the input circRNA coordinates belong to Hg38")
            hg38_matches = helpers.query_parquet(
                "test_data/annotations.parquet", "hg38", "in", coordinates, ["hg38"]
            ).to_pandas()
            hg38_matches = hg38_matches["hg38"].to_list()
            hg38_matches = sorted(list(set(hg38_matches)))

            coordinates = sorted(list(set(coordinates)))

            # Best case scenario, inputs are sound:
            if all(item in hg19_matches + hg38_matches for item in coordinates):
                # if a hg19/38 coordinate matches a key in strand_info, add the strand back to the match
                # only invoked when all coordinates are valid - nice inputs by user.
                for i, coord in enumerate(hg19_matches):
                    if coord in strand_info:
                        hg19_matches[i] = coord + ":" + strand_info[coord]
                for i, coord in enumerate(hg38_matches):
                    if coord in strand_info:
                        hg38_matches[i] = coord + ":" + strand_info[coord]

                # populate output dict.
                if hg19_matches is not None:
                    logger.debug("Input circRNA coordinates belong to Hg19")
                    database_maps["hg19"] = hg19_matches
                    del database_maps["no database (coordinates)"]
                elif hg38_matches is not None:
                    logger.debug("Input circRNA coordinates belong to Hg38")
                    database_maps["hg38"] = hg38_matches
                    del database_maps["no database (coordinates)"]

            # Fix coordinates and attempt to map for user:
            elif not all(item in hg19_matches + hg38_matches for item in coordinates):
                # Stage the database map here, then extract problem coordinates
                database_maps["hg19"] = hg19_matches
                database_maps["hg38"] = hg38_matches
                del database_maps["no database (coordinates)"]

                # Extract problem coordinates
                problem_coordinates = [
                    item
                    for item in coordinates
                    if item not in hg19_matches + hg38_matches
                ]
                logger.debug(
                    str(len(problem_coordinates))
                    + " input circRNA coordinates do not match Hg19 or Hg38. Trying 0-based and 1-based coordinates for your input."
                )

                # Add strand to OK coordinates - must be after isolating problem or else they are treated as problem
                for i, coord in enumerate(hg19_matches):
                    if coord in strand_info:
                        hg19_matches[i] = coord + ":" + strand_info[coord]

                for i, coord in enumerate(hg38_matches):
                    if coord in strand_info:
                        hg38_matches[i] = coord + ":" + strand_info[coord]

                hg19_counter = 0
                hg38_counter = 0
                for coordinate in problem_coordinates:
                    zero_based = adjust_coordinates([coordinate], "subtract")
                    one_based = adjust_coordinates([coordinate], "add")

                    zero_19 = (
                        helpers.query_parquet(
                            "test_data/annotations.parquet",
                            "hg19",
                            "in",
                            zero_based,
                            ["hg19"],
                        )
                        .to_pandas()["hg19"]
                        .drop_duplicates()
                        .tolist()
                    )
                    one_19 = (
                        helpers.query_parquet(
                            "test_data/annotations.parquet",
                            "hg19",
                            "in",
                            one_based,
                            ["hg19"],
                        )
                        .to_pandas()["hg19"]
                        .drop_duplicates()
                        .tolist()
                    )

                    if len(zero_19 + one_19) > 0:
                        hg19_counter += 1

                    zero_38 = (
                        helpers.query_parquet(
                            "test_data/annotations.parquet",
                            "hg38",
                            "in",
                            zero_based,
                            ["hg38"],
                        )
                        .to_pandas()["hg38"]
                        .drop_duplicates()
                        .tolist()
                    )
                    one_38 = (
                        helpers.query_parquet(
                            "test_data/annotations.parquet",
                            "hg38",
                            "in",
                            one_based,
                            ["hg38"],
                        )
                        .to_pandas()["hg38"]
                        .drop_duplicates()
                        .tolist()
                    )

                    if len(zero_38 + one_38) > 0:
                        hg38_counter += 1

                    # Corrected cooridinate in dict, append to database maps.
                    # ! add strand back to the corrected coordinate
                    var_dict = {"hg19": zero_19 + one_19, "hg38": zero_38 + one_38}
                    if len(var_dict["hg19"]) > 0:
                        if coordinate in strand_info:
                            var_dict["hg19"][0] = (
                                var_dict["hg19"][0] + ":" + strand_info[coordinate]
                            )
                        database_maps["hg19"].append(var_dict["hg19"][0])
                    else:
                        if coordinate in strand_info:
                            var_dict["hg38"][0] = (
                                var_dict["hg38"][0] + ":" + strand_info[coordinate]
                            )
                        database_maps["hg38"].append(var_dict["hg38"][0])

                    # Isolate corrected coordinate, use to re-write user input file
                    # {'old_coord': ['new_coord']}
                    corrected_coordinate = [list for list in var_dict.values() if list]
                    replacements[coordinate] = corrected_coordinate[0][0]

                if hg19_counter > 0 and hg38_counter == 0:
                    logger.debug("circRNA coordinates corrected and mapped to Hg19")
                elif hg19_counter == 0 and hg38_counter > 0:
                    logger.debug("circRNA coordinates corrected and mapped to Hg38")
                else:
                    logger.debug(
                        "circRNA coordinates corrected and mapped to both Hg19 and Hg38"
                    )

                # Now perform "surgery" on the input file to inject the correct coordinates for the user
                # must remain within this scope (arraystar probes wont use replacement variable)
                for row in valid_rows:
                    # key = orig, value = corrected
                    for key, values in replacements.items():
                        if row[identifier_colname] == key:
                            row[identifier_colname] = values

        # Basic info about inputs
        logger.info(
            str(len(valid_rows))
            + " circRNAs match the pycircdb database"
            )
        if len(invalid_rows) > 0:
            logger.info(
                str(len(invalid_rows))
                + " circRNA records were not found in the pycircdb database"
                )

        # Write sep files for the user. Should toggle this on/off via param
        valid_out = Path(config.outdir + "/corrected_circrna" + ext)
        invalid_out = Path(config.outdir + "/invalid_circrna" + ext)

        if len(invalid_rows) > 0:
            with invalid_out.open(mode="w", newline="") as out_handle:
                writer = csv.DictWriter(out_handle, header, delimiter=delim)
                writer.writeheader()
                for row in invalid_rows:
                    writer.writerow(row)

            logger.info(
                "Writing invalid circRNA records to "
                + str(invalid_out)
            )

        if len(valid_rows) > 0 and len(replacements) > 0:
            with valid_out.open(mode="w", newline="") as out_handle:
                writer = csv.DictWriter(out_handle, header, delimiter=delim)
                writer.writeheader()
                for row in valid_rows:
                    writer.writerow(row)

            logger.info(
                "pycircdb has fixed "
                + str(len(replacements))
                + " circRNA records and written to "
                + str(valid_out)
            )

        # Return the database map. Perform 'surgery' on it upstream within the coordinate scope
        # Remove a key if its empty
        database_maps = dict((k, v) for k, v in database_maps.items() if v)

        # output hg38 : original input so we can access downstream DB using
        # hg38, and write results using original input.
        # database maps contains all valid circRNA identifiers..

        result = defaultdict(list)

        for k, v in database_maps.items():
            if k == "hg38":
                for i in v:
                    coords = i.split(":")
                    result[coords[0:2]].append(i)
            elif k =="hg19":
                for i in v:
                    coords = i.split(":")
                    query = helpers.query_parquet("test_data/annotations.parquet", k, "in", [":".join(coords[0:2])], ["hg38", k]).to_pandas().drop_duplicates()
                    result[query["hg38"].iloc[0]].append(i)
            else:
                query = helpers.query_parquet("test_data/annotations.parquet", k, "in", v, ["hg38", k]).to_pandas().drop_duplicates()
                for row in query.iterrows():
                    result[row[1]["hg38"]].append(row[1][k])

        # convert to traditional dict:
        res = {k: v[0] for k, v in result.items()}
        return res

def ingest_others(file_in, workers, type):

    file_path = Path(file_in)
    _, ext = os.path.splitext(file_path)

    with file_path.open(newline="") as in_handle:
    
        file_dialect = sniff_format(in_handle)
        file_reader = csv.DictReader(in_handle, dialect=file_dialect)
        file_delimiter = file_dialect.delimiter
        file_header = list(file_reader.fieldnames)

        identifier_colname = file_reader.fieldnames[0]

        file_rows = []
        input_identifiers = []

        for row in file_reader:

            file_rows.append(row)
            input_identifiers.append(row[identifier_colname])

        logger.info(f'User provided {len(input_identifiers)} {type} records')

        if type == "miRNA":
            input_identifiers = [re.sub('^hsa-', '', identifier) for identifier in input_identifiers]

        database_file = match_type_to_file(type)

        matching_identifiers = helpers.check_identifiers(database_file, input_identifiers)

        identifiers_in_database = [x for x in input_identifiers if x in matching_identifiers]
        identifiers_not_in_database = [x for x in input_identifiers if x not in matching_identifiers]

        if len(identifiers_in_database) == 0:
            logger.critical(f"No valid {type} records found in the input file {file_path}'")
            return {"config": None, "report":None, "sys_exit_code": 1}

        logger.info(f'{len(identifiers_in_database)} {type}s match the pycircdb database')

        invalid_input_rows = return_invalid_input_rows(file_rows, identifiers_not_in_database, workers)

        report_invalid_input_rows = Path(config.outdir + f'/invalid_{type}_records' + ext)
        if len(invalid_input_rows) > 0:
            logger.info(f"{len(invalid_input_rows)} {type}s were not found in the pycircdb database")
            logger.info(f"Writing invalid {type} records to {report_invalid_input_rows}")
            with report_invalid_input_rows.open("w", newline="") as out_handle:
                writer = csv.DictWriter(out_handle, fieldnames=file_header, delimiter=file_delimiter)
                writer.writeheader()
                for row in invalid_input_rows:
                    writer.writerow(row)

        return matching_identifiers

databases = [
    "ArrayStar",
    "circBank",
    "circBase",
    "circAtlas",
    "no database (coordinates)",
]

def match_type_to_file(type):
    if type == "miRNA":
        file = Path("/home/barry/Desktop/pycircdb/test_data/mirna_identifiers.txt")
    elif type == "gene":
        file = Path("/home/barry/Desktop/pycircdb/test_data/mrna_identifiers.txt")
    elif type == "RBP":
        file = Path("/home/barry/Desktop/pycircdb/test_data/rbp_identifiers.txt")

    return file

def return_invalid_input_rows(file_rows, identifiers_not_in_database, workers):
        if workers > 1:
            chunk_size = len(file_rows) // mp.cpu_count()
            chunks = [file_rows[i:i + chunk_size] for i in range(0, len(file_rows), chunk_size)]
            with mp.Pool(processes=workers) as pool:
                results = pool.starmap(process_chunk, [(chunk, identifiers_not_in_database) for chunk in chunks])

            bad_entries = [item for sublist in results for item in sublist]

        else:
            bad_entries = [d for d in file_rows if any(re.search(f'{item}$', v) for item in identifiers_not_in_database for v in d.values())]

        return bad_entries

def process_chunk(chunk, list):
    return [d for d in chunk if any(item in v for item in list for v in d.values())]

def adjust_coordinates(coord_list, operation):
    adjusted_coords = []
    for coord in coord_list:
        chr_part, num_part = coord.split(":")
        start, end = num_part.split("-")
        if operation == "add":
            adjusted_start = str(int(start) + 1)
        elif operation == "subtract":
            adjusted_start = str(int(start) - 1)
        else:
            raise ValueError("Invalid operation. Choose 'add' or 'subtract'.")
        adjusted_coord = ":".join([chr_part, "-".join([adjusted_start, end])])
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
        logger.critical(f"The input file does not have a header")
        sys.exit(1)
    dialect = sniffer.sniff(peek)
    if dialect.delimiter not in [",", "\t"]:
        logger.debug(f"Single column file provided by user. Default to tab delimiter")
        dialect.delimiter = "\t"

    return dialect

if __name__ == "__main__":
    stage_inputs()
