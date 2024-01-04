import sys
import os
import rich_click as click
import pandas as pd
import pyarrow.parquet as pq
import subprocess

from pycircdb.util.helpers import query_parquet

# Set up logging
start_execution_time = time.time()
logger = config.logger

@click.group(chain=True, invoke_without_command=False)
def cli():
    """ pycircdb: A command line tool for annotating circRNAs
        using publicly available circRNA repsitories.
    """
    pass

@cli.command('interactive', short_help='Launch interactive mode')
@click.option('--launch/--no-launch', default=False, help='Toggle launching browser session.')
def interactive(launch):
    """ View database queries in interactive browser. """
    if launch:
        subprocess.run(['streamlit', 'run', 'app.py'])
    else:
        click.echo('Not launching interactive mode')

@cli.command('convert')
@click.option('-i', '--input', required=True, type=click.Path('readable'))
def convert(input):
    df = pd.read_csv(input, sep='\t', header=0, index_col=0)
    nom = sniff_nomenclature(df)
    convert_ids(nom)


def convert_ids(nom):
    #print(nom.head(4))
    parquet_column = nom.name
    parquet_to_pandas(parquet_column, nom)




if __name__ == '__main__':
    cli()