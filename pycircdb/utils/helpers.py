import os
import pyarrow
import pyarrow.parquet as pq
import dask.dataframe as dd
import dask
from dask.distributed import Client
import dask.array as da
import pandas as pd
import numpy as np
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

def query_parquet(
        parquet_file, 
        key, 
        operator, 
        value,
        columns=None) -> pyarrow.Table:
        """
        Queries parquet files

        Parameters
        ----------
        parquet_file : str
            The parquet file path.
        key : str
            The schema name to query value(s) against
        operator: str
            Conditional operator for query
        value: str, set, list, tuple
            Values to subset parquet file by
        columns: list
            Columns from parquet table to return

        Returns
        -------
        pyarrow.Table
            Pyarrow table containing subset data

        """
        # Check parquet file exists
        assert os.path.isfile(parquet_file), logger.error(f'Parquet file {parquet_file} does not exist')

        # Check the operator is valid
        valid_operators = ['=', '==', '!=', '<', '>', '<=', '>=', 'in', 'not in']
        assert operator in valid_operators, logger.error(f'{operator} is not a valid operator. Valid operators are {valid_operators}')

        # Check the value is a set, list or tuple when..
        if operator in ['in', 'not in']:
            assert isinstance(value, (set, list, tuple)), logger.error(f'value is not a set, list or tuple')

        # Check the key is valid in schema 
        assert key in pq.read_schema(parquet_file).names, logger.error(f'{key} is not a valid key in {parquet_file} schema')

        # Check the columns are valid in schema
        if columns is not None:
            for item in columns:
                assert item in pq.read_schema(parquet_file).names, logger.error(f'{item} is not a valid column in {parquet_file} schema')

        # Construct the filter using  DNF notation
        selection = [(key, operator, value)]

        # Query parquet file
        if columns:
            parquet_table = pq.read_table(parquet_file, filters=selection, columns=columns)
        else:
            parquet_table = pq.read_table(parquet_file, filters=selection)

        # debugging print statement
        #print(parquet_table.to_pandas())

        return parquet_table

def dask_query_parquet(
        parquet_file, 
        key, 
        operator, 
        value,
        columns=None) -> pyarrow.Table:
        """
        Queries parquet files

        Parameters
        ----------
        parquet_file : str
            The parquet file path.
        key : str
            The schema name to query value(s) against
        operator: str
            Conditional operator for query
        value: str, set, list, tuple
            Values to subset parquet file by
        columns: list
            Columns from parquet table to return

        Returns
        -------
        pyarrow.Table
            Pyarrow table containing subset data

        """
        # Check parquet file exists
        #assert os.path.isfile(parquet_file), logger.error(f'Parquet file {parquet_file} does not exist')

        # Check the operator is valid
        valid_operators = ['=', '==', '!=', '<', '>', '<=', '>=', 'in', 'not in']
        assert operator in valid_operators, logger.error(f'{operator} is not a valid operator. Valid operators are {valid_operators}')

        # Check the value is a set, list or tuple when..
        if operator in ['in', 'not in']:
            assert isinstance(value, (set, list, tuple)), logger.error(f'value is not a set, list or tuple')

        # Check the key is valid in schema 
        #assert key in pq.read_schema(parquet_file).names, logger.error(f'{key} is not a valid key in {parquet_file} schema')

        # Check the columns are valid in schema
        #if columns is not None:
        #    for item in columns:
        #        assert item in pq.read_schema(parquet_file).names, logger.error(f'{item} is not a valid column in {parquet_file} schema')

        # Construct the filter using  DNF notation
        selection = [(key, operator, value)]

        # Query parquet file
        if columns:
            parquet_table = dd.read_parquet(parquet_file, filters=selection, columns=columns, split_row_groups="adaptive")
        else:
            parquet_table = dd.read_parquet(parquet_file, filters=selection, split_row_groups="adaptive")

        # debugging print statement
        #print(parquet_table.to_pandas())

        return parquet_table


def process_db(df, identifier):
    mask = df['ID'].isin(identifier)
    subset = df[mask].dropna()
    return subset['ID'].tolist()


def check_identifiers(
        id_file,
        identifiers=None
        ) -> list:
        """
        miRNA, RBP and HGNC symbols are in a one column txt file. 
        Checking if the users inputs match the databses can be done
        using these smaller text files. 
        """

        df = pd.read_csv(id_file, names=["ID"], dtype=str)
        results = []
        x = process_db(df, identifiers)
        results.extend(x)

        return results

