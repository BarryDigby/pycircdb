import argparse
import pyarrow.parquet as pq

# Create the parser
parser = argparse.ArgumentParser(description='Read a parquet file.')

# Add the arguments
parser.add_argument('FilePath', metavar='filepath', type=str, help='the path to the parquet file')

# Parse the arguments
args = parser.parse_args()

# Read the parquet file
parquet_file = pq.ParquetFile(args.FilePath)
print("\nSchema using pyarrow:")
print(parquet_file.schema)

