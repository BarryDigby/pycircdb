import sys
import os
import rich
import rich_click as click
import pandas as pd
import pyarrow.parquet as pq
import subprocess
import logging
import time

from .utils.helpers import query_parquet
from .utils.ingest import check_circrna
from .utils.annotate import annotate_circRNAs
from .utils import config, log, util_functions, report

# Set up logging
start_execution_time = time.time()
logger = config.logger

""" 
    From ~/Desktop/pycircdb
    Run: python3.10 -m pip install -e . 
    pycircdb --help should now work
"""
# Configuration for rich-click CLI help
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True
click.rich_click.HEADER_TEXT = f"[dark_orange]>>>[/] [bold][link=https://github.com/BarryDigby/pycircdb]pycircdb[/link][/] [dim]| v{config.version}"
click.rich_click.FOOTER_TEXT = "See [link=https://github.com/BarryDigby/pycircdb]https://github.com/BarryDigby/pycircdb[/] for more details."
click.rich_click.ERRORS_SUGGESTION = f"This is pycircdb [cyan]v{config.version}[/]\nFor more help, run '[yellow]pycircdb --help[/]' or visit [link=https://github.com/BarryDigby/pycircdb]https://github.com/BarryDigby/pycircdb[/]"
click.rich_click.STYLE_ERRORS_SUGGESTION = ""
click.rich_click.OPTION_GROUPS = {
    "pycircdb": [
        {
            "name": "Main Options",
            "options": ["--circRNA", "--miRNA", "--mRNA", "--RBP", "--interactive"],
            "help": "here is help",
        },
        {
            "name": "Options",
            "options": [
                "--annotate",
                "--network",
                "--ceRNA_network",
                "--outdir",
                "--verbose",
                "--quiet",
            ],
        },
        {
            "name": "Filtering Options",
            "options": ["--circRNA_algorithm", "--circRNA_miRNA_DB", "--miRNA_mRNA_DB"],
        },
    ]
}


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "--circRNA",
    "circRNA",
    type=click.Path(exists=True, readable=True),
    help="Input circRNA file conatining circRNA identifiers in the first column (i.e rownames).\n\nAccepted circRNA identifiers: Hg19, Hg38, ArrayStar, CircBank, CircBase, circAtlas",
)
@click.option(
    "--miRNA",
    "miRNA",
    type=click.Path(exists=True, readable=True),
    help="Input miRNA file containing miRNA identfiers in the first column (i.e rownames).\n\nAccepted miRNA identifiers: latest miRBase release (hsa-miR-106a-5p/miR-106a-5p)",
)
@click.option(
    "--mRNA",
    "mRNA",
    type=click.Path(exists=True, readable=True),
    help="Input gene file containing gene identifiers in the first column (i.e rownames).\n\nAccepted gene identifiers: HGNC symbols (REG4, PTEN, etc)",
)
@click.option(
    "--RBP",
    "RBP",
    type=click.Path(exists=True, readable=True),
    help="Input RNA binding protein file containing RBP identifiers in the first column (i.e rownames).\n\nAccepted RBP identifiers: Gene symbols (ELAVL1, AGO2, etc)",
)
@click.option(
    "--interactive",
    "interactive",
    is_flag=True,
    help="Deploy interactive Streamlit HTML to view results",
)
@click.option(
    "--annotate",
    "annotate",
    is_flag=True,
    help="Annotate circRNAs passed via `--circRNA` parameter.",
)
@click.option(
    "--network", "network", is_flag=True, help="Construct circRNA-miRNA-mRNA network"
)
@click.option(
    "--ceRNA-network",
    "ceRNA_network",
    is_flag=True,
    help="Construct circRNA-miRNA-mRNA ceRNA network",
)
@click.option("--outdir", "outdir", type=str, help="Output directory")
@click.option(
    "--circRNA-algorithm",
    "circRNA_algorithm",
    type=click.Choice(["CIRCexplorer2", "circRNA_finder", "CIRI", "find_circ"]),
    multiple=True,
    help="circRNA detection algorithm",
)
@click.option(
    "--circRNA-miRNA-DB",
    "circRNA_miRNA_DB",
    type=click.Choice(["miRanda", "TargetScan"]),
    multiple=True,
    help="circRNA-miRNA database",
)
@click.option(
    "--miRNA-mRNA-DB",
    "miRNA_mRNA_DB",
    type=click.Choice(["miRDB", "miRTarBase", "TarBase", "TargetScan"]),
    multiple=True,
    help="miRNA-mRNA database",
)
@click.option("--data-dir", "data_dir", is_flag=True, help="Data directory")
@click.option("--verbose", "verbose", count=True, default=0, help="Verbose output")
@click.option("--quiet", "quiet", is_flag=True, help="Only show Log warnings")
@click.version_option(config.version, prog_name="pycircdb")
def run_cli(**kwargs):
    """
    pycircdb: a command line tool that pools publicly available circRNA databases. . . .

    At a minimum, the package accepts as input a list of circRNAs (one per line) and returns a dataframe of annotations for the user.

    '[blue bold]$ pycircdb --input circrna.txt[/]'
    """
    # Pass on to a regular function that can be used easily without click
    pycircdb_run = run(**kwargs)

    # End execution using the exit code returned from pycircdb
    sys.exit(pycircdb_run["sys_exit_code"])


# Main function that runs pycircdb. Available to use within an interactive Python environment
def run(
    circRNA=None,
    miRNA=None,
    mRNA=None,
    RBP=None,
    annotate=True,
    network=False,
    ceRNA_network=False,
    outdir=None,
    circRNA_algorithm=None,
    circRNA_miRNA_DB=None,
    miRNA_mRNA_DB=None,
    verbose=0,
    no_ansi=False,
    quiet=False,
    **kwargs,
):
    """
    Main run function for pycircdb
    """

    # Set up logging level
    loglevel = log.LEVELS.get(min(verbose, 1), "INFO")
    if quiet:
        loglevel = "WARNING"
        config.quiet = True
    log.init_log(logger, loglevel=loglevel, no_ansi=no_ansi)

    # Set up rich console
    console = rich.console.Console(
        stderr=True,
        highlight=False,
        force_terminal=util_functions.force_term_colors(),
        color_system=None if no_ansi else "auto",
    )
    console.print(
        f"\n[dark_orange]>>>[/] [bold][link=https://github.com/BarryDigby/pycircdb]pycircdb[/link][/] [dim]| v{config.version}\n"
    )

    # Versions, commands used etc for debug
    logger.debug(f"This is pycircdb v{config.version}")
    logger.debug("Running Python " + sys.version.replace("\n", " "))
    # Log the command used to launch pycircdb
    report.init()
    report.pycircdb_command = " ".join(sys.argv)
    logger.debug(f"Command used: {report.pycircdb_command}")

    # Set up key variables (overwrite config vars from command line)
    if outdir is not None:
        config.output_dir = outdir

    # Delete and use config going forward
    del outdir

    # Set up results directory
    if not os.path.exists(config.output_dir):
        os.makedirs(config.output_dir)

    # Assess user inputs
    circrnas = check_circrna(circRNA) # returns dict of {'database': ['identifier'], 'database2': ['identifier]}
    
    # does the user want to annotate the circRNAs?
    if annotate:
        logger.info("Annotating circRNAs")
        annotate_circRNAs(circrnas)
    else:
        logger.info("Skipping circRNA annotation")

    sys_exit_code = 0
    # return dict for pycircdb_run
    return {"sys_exit_code": sys_exit_code}
