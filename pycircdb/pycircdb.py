import sys
import os
import rich
import rich_click as click
import pandas as pd
import pyarrow.parquet as pq
import subprocess
import logging
import time
from itertools import permutations
from .utils.parameter_check import ingest
from .utils.ingest import stage_inputs
from .utils.filter_circrnas import filter_circrna
from .utils.annotate import annotate_circrna
from .utils.filter_mirnas import filter_mirna
from .utils.filter_genes import filter_gene
from .utils.build_network import initialise_network


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
            "name": "Input options",
            "options": ["--circrna", "--mirna", "--gene", "--rbp"],
        },
        {
            "name": "Modules",
            "options": [
                "--annotate-circrnas",
                "--mirna-targets",
                "--rna-binding-proteins",
                "--cerna-network",
            ],
        },
        {
            "name": "Filtering options",
            "options": [
                "--circrna-algorithm",
                "--circrna-set-logic",
                "--mirna-algorithm",
                "--mirna-set-logic",
                "--mirna-type",
                "--mirna-mfe",
                "--mirna-score",
                "--gene-database",
                "--gene-set-logic",
            ],
        },
        {
            "name": "Options",
            "options": ["--workers", "--outdir", "--verbose", "--quiet", "--version", "--help"],
        },
    ]
}

def default_none(ctx, param, value):
    if len(value) == 0:
        return None
    else:
        return value

@click.command(context_settings=dict(help_option_names=["-h", "--help"]))

# Input options:

@click.option(
    "--circrna",
    "circrna_file",
    type=click.Path(exists=True, readable=True),
    help="Input circRNA file conatining circRNA identifiers in the first column.\n\nAccepted circRNA identifiers: Hg19, Hg38, ArrayStar, CircBank, CircBase, circAtlas",
)
@click.option(
    "--mirna",
    "mirna_file",
    type=click.Path(exists=True, readable=True),
    help="Input miRNA file containing miRNA identfiers in the first column.\n\nAccepted miRNA identifiers: latest miRBase release (hsa-miR-106a-5p, miR-106a-5p)",
)
@click.option(
    "--gene",
    "gene_file",
    type=click.Path(exists=True, readable=True),
    help="Input gene file containing gene identifiers in the first column.\n\nAccepted gene identifiers: HGNC symbols (REG4, PTEN, etc)",
)
@click.option(
    "--rbp",
    "rbp_file",
    type=click.Path(exists=True, readable=True),
    help="Input RNA binding protein file containing RBP identifiers in the first column.\n\nAccepted RBP identifiers: Gene symbols (ELAVL1, AGO2, etc)",
)

# pycircdb modules:

@click.option(
    "--annotate-circrnas",
    "annotate_circrnas",
    is_flag=True,
    help="Annotate circRNAs passed via `--circRNA` parameter.",
)
@click.option(
    "--mirna-targets",
    "mirna_targets",
    is_flag=True,
    help="Construct circRNA-miRNA network. Requires `--circRNA` or `--miRNA` file. Both files can be supplied to initialise the network.",
)
@click.option(
    "--rna-binding-proteins",
    "rna_binding_proteins",
    is_flag=True,
    help="Construct circRNA-RBP network. Requires `--circRNA` or `--RBP` file. Both files can be supplied to initialise the network.",
)
@click.option(
    "--cerna-network",
    "cerna_network",
    is_flag=True,
    help="Construct circRNA-miRNA-mRNA network. Requires `--circRNA`, `--miRNA` or `--mRNA` file.\n\nAny combination of files can be supplied to initialise the network.",
)

# Network filtering options:


@click.option(
    "--circrna-algorithm",
    "-ca",
    "circrna_algorithm",
    type=click.Choice(config.circrna_algorithm),
    multiple=True,
    callback=default_none,
    help="Subset network using circRNAs detected by circexplorer2, circrna_finder, ciri or find_circ.\n\nCan be specified multiple times: --circRNA-algorithm ciri --ca circexplorer2",
)
@click.option(
    "--circrna-set-logic",
    "circrna_set_logic",
    type=click.Choice(["AND", "OR"]),
    default="AND",
    help="Set logic to apply when multiple circRNA detection algorithms selected using `--circrna-algorithm`.\n\nDefault: AND.",
)
@click.option(
    "--mirna-algorithm",
    "-ma",
    "mirna_algorithm",
    type=click.Choice(config.mirna_algorithm),
    multiple=True,
    callback=default_none,
    help="Subset network using circRNA-miRNAs interactions predicted by miRanda and/or TargetScan.",
)
@click.option(
    "--mirna-set-logic",
    "mirna_set_logic",
    type=click.Choice(["AND", "OR"]),
    default="AND",
    help="Set logic to apply when multiple miRNA detection algorithms selected using `--mirna-algorithm`.\n\nDefault: AND.",
)
@click.option(
    "--mirna-type",
    "-mt",
    "mirna_type",
    type=click.Choice(config.mirna_type),
    multiple=True,
    callback=default_none,
    help="Subset network using predicted circRNA-miRNA TargetScan MRE site types.",
)
@click.option(
    "--mirna-mfe",
    "-mfe",
    "mirna_mfe",
    type=click.FloatRange(min=-62.0, max=-0.41),
    help="Subset network using predicted circRNA-miRNA interaction minimum free energy (MFE) threshold.\n\nRange: -62.0, -0.41\n\nPlease note this only applies to circRNA-miRNA interactions predicted by miRanda.",
)
@click.option(
    "--mirna-score",
    "-ms",
    "mirna_score",
    type=click.FloatRange(min=140.0, max=220.0),
    help="Subset network circRNA-miRNA miRanda interaction scores.\n\nRange: 140.0, 220.0\n\nPlease note this only applies to circRNA-miRNA interactions predicted by miRanda.",
)
@click.option(
    "--gene-database",
    "-gdb",
    "gene_database",
    type=click.Choice(config.gene_database),
    multiple=True,
    callback=default_none,
    help="Subset network using miRNA-mRNA interactions predicted by databases.\n\nValidated miRNA-mRNA interaction databases: miRecords, miRTarBase, TarBase.\n\nPredicted miRNA-mRNA interaction databases: DIANA, ElMMo, MicroCosm, miRanda, miRDB, PicTar, PITA TargetScan",
)
@click.option(
    "--gene-set-logic",
    "gene_set_logic",
    type=click.Choice(["AND", "OR"]),
    default="AND",
    help="Set logic to apply when multiple gene databases selected using `--gene-database`.\n\nDefault: AND.",
)

# pycircdb options:

@click.option("--workers", "-c", "workers", type=int, help="Number of cores")
@click.option("--outdir", "-o", "outdir", required=True, type=str, help="Output directory")
@click.option("--verbose", "verbose", count=True, default=0, help="Verbose output")
@click.option("--quiet", "quiet", is_flag=True, help="Only show Log warnings")
@click.option("--version", "-v", "version", is_flag=True, help="Show version and exit")
@click.version_option(config.version, prog_name="pycircdb")

def run_cli(**kwargs):
    """
    pycircdb: a command line tool that pools publicly available circRNA databases. . . .

    At a minimum, the package accepts as input a list of circRNAs (one per line) and returns a dataframe of annotations for the user.

    '[blue bold]$ pycircdb --circRNA circrna.txt[/]'
    """
    # Pass on to a regular function that can be used easily without click
    pycircdb_run = run(**kwargs)

    # End execution using the exit code returned from pycircdb
    sys.exit(pycircdb_run["sys_exit_code"])


# Main function that runs pycircdb. Available to use within an interactive Python environment
def run(
    circrna_file=None,
    mirna_file=None,
    gene_file=None,
    rbp_file=None,
    annotate_circrnas=False,
    mirna_targets=False,
    rna_binding_proteins=False,
    cerna_network=False,
    circrna_algorithm=None,
    circrna_set_logic=None,
    mirna_algorithm=None,
    mirna_set_logic=None,
    mirna_type=None,
    mirna_mfe=None,
    mirna_score=None,
    gene_database=None,
    gene_set_logic=None,
    workers=None,
    outdir=None,
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
    # if not $var logic is used below to set undefined variables to None (otherwise populated by config defaults)
    if outdir is not None:
        config.outdir = outdir
    if workers is not None:
        config.workers = workers
    if annotate_circrnas is not None:
        config.annotate_circrnas = annotate_circrnas
    if not circrna_algorithm:
        config.circrna_algorithm = circrna_algorithm
    else:
        config.circrna_algorithm = circrna_algorithm
    if circrna_set_logic is not None:
        config.circrna_set_logic = circrna_set_logic
    if not mirna_algorithm:
        config.mirna_algorithm = mirna_algorithm
    else:
        config.mirna_algorithm = mirna_algorithm
    if mirna_set_logic is not None:
        config.mirna_set_logic = mirna_set_logic
    if not mirna_type:
        config.mirna_type = mirna_type
    else:
        config.mirna_type = mirna_type
    if mirna_mfe is not None:
        config.mirna_mfe = mirna_mfe
    if mirna_score is not None:
        config.mirna_score = mirna_score
    if not gene_database:
        config.gene_database = gene_database
    else:
        config.gene_database = gene_database
    if gene_set_logic is not None:
        config.gene_set_logic = gene_set_logic
    


    # Delete and use config going forward
    del outdir
    del workers
    del annotate_circrnas
    del circrna_algorithm
    del circrna_set_logic
    del mirna_algorithm
    del mirna_set_logic
    del mirna_type
    del mirna_mfe
    del mirna_score
    del gene_database
    del gene_set_logic

    # Set up results directory
    if not os.path.exists(config.outdir):
        os.makedirs(config.outdir)

    # Assess user inputs
    circrna_dict, mirna_list, gene_list, rbp_list = stage_inputs(circrna_file, mirna_file, gene_file, rbp_file, config.workers)

    # Filter biotypes using user inputs - broken into seperate functions to simplify execution
    annotated_df = None
    if circrna_dict is not None and config.circrna_algorithm is not None:
        filtered_circrnas, annotated_df = filter_circrna(circrna_dict, config.circrna_algorithm, config.circrna_set_logic)
    else:
        filtered_circrnas = circrna_dict
    
    # Return annotated circRNAs if user wants them 
    if config.annotate_circrnas:
        annotate_circrna(filtered_circrnas, annotated_df)

    # miRNA filtering.
    if mirna_list is not None:
        filtered_mirnas = filter_mirna(mirna_list, config.mirna_algorithm, config.mirna_set_logic, config.mirna_type, config.mirna_mfe, config.mirna_score, config.workers)

    # gene filtering
    # basic triggers for gene filtering
    if gene_list is not None:
        filtered_genes = filter_gene(gene_list, config.gene_database, config.gene_set_logic, config.workers)

    # Build network files
    #foo = initialise_network(filtered_circrnas, mirna_list, gene_list, rbp_list, config.workers)
        
    # dont delete me
    sys_exit_code = 0
    # return dict for pycircdb_run
    return {"sys_exit_code": sys_exit_code}

