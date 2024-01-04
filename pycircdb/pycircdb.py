import sys
import os
import rich_click as click
import pandas as pd
import pyarrow.parquet as pq
import subprocess

from .utils.helpers import query_parquet
from .utils import config

""" 
    From ~/Desktop/pycircdb
    Run: python3.10 -m pip install -e . 
    pycircdb --help should now work
"""
# Configuration for rich-click CLI help
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True
click.rich_click.HEADER_TEXT = (
    f"[dark_orange]///[/] [bold][link=https://github.com/BarryDigby/pycircdb]pycircdb[/link][/] :mag: :loop: :curly_loop: [dim]| v{config.version}"
)
click.rich_click.FOOTER_TEXT = "See [link=https://github.com/BarryDigby/pycircdb]https://github.com/BarryDigby/pycircdb[/] for more details."
click.rich_click.ERRORS_SUGGESTION = f"This is pycircdb [cyan]v{config.version}[/]\nFor more help, run '[yellow]pycircdb --help[/]' or visit [link=https://github.com/BarryDigby/pycircdb]https://github.com/BarryDigby/pycircdb[/]"
click.rich_click.STYLE_ERRORS_SUGGESTION = ""
click.rich_click.OPTION_GROUPS = {
    "pycircdb": [
        {
            "name": "Main options",
            "options": [
                "--input",
                "--reference"
            ]
        }
    ]
}

@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option("-i", "--input", type=click.Path(exists=True, readable=True), required=True, help="Input circRNA list.")
@click.option("-r", "--reference", type=click.Choice(["hg19", "hg38"]), default="hg19", required=True, help="Reference genome build.")
@click.option("-o", "--outdir", type=str, help="Output results to the specified output directory.")
@click.option("-n", "--filename", type=str, help="Append prefix to output filenames.")
@click.version_option(config.version, prog_name="pycircdb")

def run_cli(**kwargs):
    """
    pycircdb: a command line tool for rich annotation of circRNAs using publicly available circRNA repsitories.

    At a minimum, the package accepts as input a list of circRNAs (one per line) and returns a dataframe of annotations for the user. 

    '[blue bold]$ pycircdb --input circrna.txt[/]'
    """
    # Pass on to a regular function that can be used easily without click
    pycircdb_run = run(**kwargs)

    # End execution using the exit code returned from pycircdb
    sys.exit(pycircdb_run["sys_exit_code"])


# Main function that runs pycircdb. Available to use within an interactive Python environment
def run(
    input=None,
    reference=None,
    outdir=None,
    filename=None,
    **kwargs):
    """
    Main run function for pycircdb
    """

    print(f'{input}, {reference}')

    sys_exit_code = 0

    # return dict for pycircdb_run
    return {"sys_exit_code": sys_exit_code}
