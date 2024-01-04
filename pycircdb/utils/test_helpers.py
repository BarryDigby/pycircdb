import os
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from pycircdb.util.helpers import query_parquet

@pytest.fixture
def parquet_file(tmpdir):
    # Create a temporary parquet file for testing
    parquet_file_path = str(tmpdir.join("test.parquet"))
    table = pa.Table.from_pandas(pd.DataFrame({"key": [1, 2, 3], "value": ["a", "b", "c"]}))
    pq.write_table(table, parquet_file_path)
    return parquet_file_path

def test_query_parquet_existing_file(parquet_file):
    # Test querying an existing parquet file
    result = query_parquet(parquet_file, "key", "=", 2)
    assert len(result) == 1
    assert result[0]["key"] == 2
    assert result[0]["value"] == "b"

def test_query_parquet_non_existing_file():
    # Test querying a non-existing parquet file
    with pytest.raises(AssertionError):
        query_parquet("non_existing.parquet", "key", "=", 2)

def test_query_parquet_invalid_operator(parquet_file):
    # Test querying with an invalid operator
    with pytest.raises(AssertionError):
        query_parquet(parquet_file, "key", "invalid_operator", 2)

def test_query_parquet_invalid_value(parquet_file):
    # Test querying with an invalid value
    with pytest.raises(AssertionError):
        query_parquet(parquet_file, "key", "in", "invalid_value")

def test_query_parquet_invalid_key(parquet_file):
    # Test querying with an invalid key
    with pytest.raises(AssertionError):
        query_parquet(parquet_file, "invalid_key", "=", 2)