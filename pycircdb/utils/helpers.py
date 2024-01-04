import os
import pyarrow
import pyarrow.parquet as pq
import pandas as pd

def query_parquet(
        parquet_file, 
        key, 
        operator, 
        value) -> pyarrow.Table:
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

        Returns
        -------
        pyarrow.Table
            Pyarrow table containing subset data

        """
        # Check parquet file exists
        assert os.path.isfile(parquet_file), f'Parquet file {parquet_file} does not exist'

        # Check the operator is valid
        valid_operators = ['=', '==', '!=', '<', '>', '<=', '>=', 'in', 'not in']
        assert operator in valid_operators, f'{operator} is not a valid operator. Valid operators are {valid_operators}'

        # Check the value is a set, list or tuple when..
        if operator in ['in', 'not in']:
            assert isinstance(value, (set, list, tuple)), f'value is not a set, list or tuple'

        # Check the key is valid in schema 
        assert key in pq.read_schema(parquet_file).names, f'{key} is not a valid key in {parquet_file} schema'

        # Construct the filter using  DNF notation
        selection = [(key, operator, value)]

        # Query parquet file
        parquet_table = pq.read_table(parquet_file, filters=selection)

        return parquet_table


def detect_circrna_(input):
    circ_id = input.index.to_series()
    prefixes = {
        'ASCRP': 'Arraystar_ID',
        'hsa_circ_': 'circBase_ID',
        'hsa-': 'circAtlas_ID',
        'hsa_circ[a-zA-Z]\w+': 'circBank_ID',
        'chr': 'Position'
    }
    
    for prefix, identifier in prefixes.items():
        if circ_id.str.startswith(prefix).any():
            if circ_id.str.startswith(prefix).all():
                print(f'All {identifier} IDs')
                circ_id.name = identifier
                return circ_id
            else:
                sys.exit('Input circRNA matrix contains a mixture of IDs as rownames')
    
    sys.exit('Input circRNA matrix contains unrecognized rownames')