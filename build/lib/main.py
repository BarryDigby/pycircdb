import json
import shutil
from pathlib import Path
import rich_click as click
from typing import Tuple, List, Optional
from config import create_config, load_config, print_config_panel

# Workflow stuff
from hamilton import driver
from hamilton.execution import executors

from utils.connect_s3.download_annotation_tables import fetch_annotation_tables
from utils.connect_s3.download_sequence_tables import fetch_sequence_tables
from utils.connect_s3.download_mirna_tables import fetch_mirna_tables
from utils.connect_s3.download_rbp_tables import fetch_rbp_tables

import utils.detect_inputs.detect_inputs_driver as instantiate_lookup_driver
import utils.annotate.annotate_driver as annotation_driver
import utils.fasta.sequence_driver as sequence_driver
import utils.mirna.mirna_driver as mirna_driver
import utils.rbp.rbp_driver as rbp_driver
from rich.console import Console

console = Console(stderr=True, highlight=False)

@click.group(chain=True, context_settings=dict(help_option_names=['-h', '--help']))
@click.option(
    "-c",
    "--config",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    required=False,
    help="Path to the JSON config file containing workflow parameters."
)
@click.option(
    "-v",
    "--verbose",
    type=click.IntRange(0, 2),
    default=1,
    help="Verbosity level: 0 (silent), 1 (high-level, default), 2 (all outputs)."
)
@click.pass_context
def cli(ctx, config, verbose):
    """Main CLI tool."""
    ctx.ensure_object(dict)
    
    ctx.obj['verbose'] = verbose
    
    if config:
        cfg = load_config(config, verbose=verbose)
        if not cfg.get('samples'):
            raise click.UsageError("Configuration file must contain a 'samples' dictionary.")
        cfg.update({'verbose': verbose})
        ctx.obj['cfg'] = cfg
    else:
        ctx.obj['cfg'] = None
    ctx.obj['lookup_dict'] = None


@cli.result_callback()
@click.pass_context
def process_pipeline(ctx, processors, config, verbose):
    """Execute all processors returned by subcommands after parsing args."""
    cfg = ctx.obj.get('cfg')
    if cfg:
        # Populate missing config for printing correctly
        if 'annotate_databases' not in cfg:
            cfg['annotate_databases'] = ['arraystar', 'circbank', 'circbase', 'circpedia', 'circrna_db', 'cscd', 'exorbase']
        if 'fasta_databases' not in cfg:
            cfg['fasta_databases'] = ['arraystar', 'circbank', 'circbase', 'circpedia', 'circrna_db', 'cscd']
        if 'mirna_algorithms' not in cfg:
            cfg['mirna_algorithms'] = ['miranda', 'pita', 'targetscan']
            
        if verbose >= 2:
            print_config_panel(cfg, config)
            
    for processor in processors:
        processor()

@cli.command('annotate')
@click.option(
    "-d",
    "--database",
    type=str,
    required=False,
    default="arraystar,circbank,circbase,circpedia,circRNA_DB,CSCD,exorbase",
    show_default=True,
    help="Comma-separated list of databases to use."
)
@click.pass_context
def annotate(ctx, database):
    """Annotate circRNAs using a JSON configuration file."""
    cfg = ctx.obj.get('cfg')
    if not cfg:
        raise click.UsageError("A config file must be provided via -c/--config before subcommands (e.g., main.py -c config.json annotate)")

    if database:
        valid_dbs = {'arraystar', 'circbank', 'circbase', 'circpedia', 'circrna_db', 'cscd', 'exorbase'}
        parsed_dbs = [d.strip().lower() for d in database.split(',')]
        invalid_dbs = [d for d in parsed_dbs if d not in valid_dbs]
        if invalid_dbs:
            raise click.BadParameter(f"Invalid databases provided: {', '.join(invalid_dbs)}. Valid options are: {', '.join(sorted(valid_dbs))}")
        cfg['annotate_databases'] = parsed_dbs
    else:
        cfg['annotate_databases'] = ["arraystar", "circbank", "circbase", "circpedia", "circrna_db", "cscd", "exorbase"]

    def processor():
        lookup_dict = ctx.obj.get('lookup_dict')
        ctx.obj['lookup_dict'] = run_annotation(lookup_dict=lookup_dict, **cfg)
    return processor


@cli.command('fasta')
@click.option(
    "-d",
    "--database",
    type=str,
    required=False,
    default="arraystar,circbank,circbase,circpedia,circRNA_DB,CSCD",
    show_default=True,
    help="Comma-separated list of databases to use."
)
@click.pass_context
def fasta(ctx, database):
    """Output circRNA sequences in FASTA format."""
    cfg = ctx.obj.get('cfg')
    if not cfg:
        raise click.UsageError("A config file must be provided via -c/--config")

    if database:
        valid_dbs = {'arraystar', 'circbank', 'circbase', 'circpedia', 'circrna_db', 'cscd'}
        parsed_dbs = [d.strip().lower() for d in database.split(',')]
        invalid_dbs = [d for d in parsed_dbs if d not in valid_dbs]
        if invalid_dbs:
            raise click.BadParameter(f"Invalid databases provided: {', '.join(invalid_dbs)}. Valid options are: {', '.join(sorted(valid_dbs))}")
        cfg['fasta_databases'] = parsed_dbs
    else:
        cfg['fasta_databases'] = ["arraystar", "circbank", "circbase", "circpedia", "circrna_db", "cscd"]

    def processor():
        lookup_dict = ctx.obj.get('lookup_dict')
        ctx.obj['lookup_dict'] = run_fasta(lookup_dict=lookup_dict, **cfg)
    return processor


@cli.command('mirna')
@click.option(
    "-a",
    "--algorithm",
    type=str,
    required=False,
    default="miRanda,PITA,TargetScan",
    show_default=True,
    help="Comma-separated list of algorithms to use."
)
@click.pass_context
def mirna(ctx, algorithm):
    """Output miRNA interactions for identified circRNAs."""
    cfg = ctx.obj.get('cfg')
    if not cfg:
        raise click.UsageError("A config file must be provided via -c/--config")

    if algorithm:
        valid_algs = {'miranda', 'pita', 'targetscan'}
        parsed_algs = [a.strip() for a in algorithm.split(',')]
        invalid_algs = [a for a in parsed_algs if a.lower() not in valid_algs]
        if invalid_algs:
            raise click.BadParameter(f"Invalid algorithms provided: {', '.join(invalid_algs)}. Valid options are: miRanda, PITA, TargetScan")
        # Keep original case for 'contains' check, or use lowercase for case-insensitive check
        cfg['mirna_algorithms'] = parsed_algs
    else:
        cfg['mirna_algorithms'] = ["miRanda", "PITA", "TargetScan"]

    def processor():
        lookup_dict = ctx.obj.get('lookup_dict')
        run_mirna(lookup_dict=lookup_dict, **cfg)
    return processor


@cli.command('rbp')
@click.pass_context
def rbp(ctx):
    """Output RBP interactions for identified circRNAs."""
    cfg = ctx.obj.get('cfg')
    if not cfg:
        raise click.UsageError("A config file must be provided via -c/--config")

    def processor():
        lookup_dict = ctx.obj.get('lookup_dict')
        run_rbp(lookup_dict=lookup_dict, **cfg)
    return processor

def run_annotation(lookup_dict=None, **kwargs):
    """Run the annotation workflow."""

    # Lookup tables
    if lookup_dict is None:
        lookup_dict = instantiate_lookup_driver.instantiate_driver(kwargs, verbose=kwargs.get('verbose', 1))

    tmp_dir = kwargs.get("global_parameters", {}).get("tmp_dir", "tmp")

    # Filter lookup dictionary if databases option is provided
    databases = kwargs.get("annotate_databases")
    if databases:
        filtered_lookup_dict = {}
        for sample, db_dict in lookup_dict.items():
            filtered_lookup_dict[sample] = {k: v for k, v in db_dict.items() if k.lower() in databases}
    else:
        filtered_lookup_dict = lookup_dict

    # Pull annotation tables
    annotation_tables = fetch_annotation_tables(filtered_lookup_dict, tmp_dir_path=tmp_dir, verbose=kwargs.get("verbose", 1))

    # Annotate + write to file in parallel
    dr = (
        driver.Builder()
        .enable_dynamic_execution(allow_experimental_mode=True)
        .with_local_executor(executors.SynchronousLocalTaskExecutor())
        .with_remote_executor(executors.MultiThreadingExecutor(max_tasks=kwargs.get("global_parameters", {}).get("max_tasks", 1)))
        .with_modules(annotation_driver)
        .build()
    )
    
    dr.execute(
        ['close_annotation'],
        inputs={'config': kwargs, 'annotation_tables': annotation_tables, 'lookup_results': filtered_lookup_dict}
    )
    
    return lookup_dict


def run_fasta(lookup_dict=None, **kwargs):
    """Generate FASTA output from circRNA sequences.
    
    Args:
        lookup_dict: Optional pre-computed lookup results. If None, will be generated from scratch.
        **kwargs: Configuration parameters.
    """
    # Generate lookup tables if not provided (i.e., fasta running standalone)
    if lookup_dict is None:
        lookup_dict = instantiate_lookup_driver.instantiate_driver(kwargs, verbose=kwargs.get("verbose", 1))

    databases = kwargs.get("fasta_databases")
    if databases:
        filtered_lookup_dict = {}
        for sample, db_dict in lookup_dict.items():
            filtered_lookup_dict[sample] = {k: v for k, v in db_dict.items() if k.lower() in databases}
    else:
        filtered_lookup_dict = lookup_dict

    tmp_dir = kwargs.get("global_parameters", {}).get("tmp_dir", "tmp")
    sequence_tables = fetch_sequence_tables(filtered_lookup_dict, tmp_dir_path=tmp_dir, verbose=kwargs.get("verbose", 1))

    dr = (
        driver.Builder()
        .enable_dynamic_execution(allow_experimental_mode=True)
        .with_local_executor(executors.SynchronousLocalTaskExecutor())
        .with_remote_executor(executors.MultiThreadingExecutor(max_tasks=kwargs.get("global_parameters", {}).get("max_tasks", 1)))
        .with_modules(sequence_driver)
        .build()
    )

    dr.execute(
        ['close_sequence'],
        inputs={'config': kwargs, 'lookup_dict': filtered_lookup_dict, 'sequence_tables': sequence_tables}
    )

    return lookup_dict


def run_mirna(lookup_dict=None, **kwargs):
    """Output miRNA interactions for identified circRNAs.

    Args:
        lookup_dict: Optional pre-computed lookup results. If None, will be generated from scratch.
        **kwargs: Configuration parameters.
    """
    if lookup_dict is None:
        lookup_dict = instantiate_lookup_driver.instantiate_driver(kwargs, verbose=kwargs.get("verbose", 1))


    tmp_dir = kwargs.get("global_parameters", {}).get("tmp_dir", "tmp")
    mirna_tables = fetch_mirna_tables(lookup_dict, tmp_dir_path=tmp_dir, verbose=kwargs.get("verbose", 1))

    dr = (
        driver.Builder()
        .enable_dynamic_execution(allow_experimental_mode=True)
        .with_local_executor(executors.SynchronousLocalTaskExecutor())
        .with_remote_executor(executors.MultiThreadingExecutor(max_tasks=kwargs.get("global_parameters", {}).get("max_tasks", 1)))
        .with_modules(mirna_driver)
        .build()
    )
    dr.execute(
        ['close_mirna'],
        inputs={'config': kwargs, 'lookup_dict': lookup_dict, 'mirna_tables': mirna_tables}
    )


def run_rbp(lookup_dict=None, **kwargs):
    """Output RBP interactions for identified circRNAs.

    Args:
        lookup_dict: Optional pre-computed lookup results. If None, will be generated from scratch.
        **kwargs: Configuration parameters.
    """
    if lookup_dict is None:
        lookup_dict = instantiate_lookup_driver.instantiate_driver(kwargs, verbose=kwargs.get("verbose", 1))


    tmp_dir = kwargs.get("global_parameters", {}).get("tmp_dir", "tmp")
    rbp_tables = fetch_rbp_tables(lookup_dict, tmp_dir_path=tmp_dir, verbose=kwargs.get("verbose", 1))

    dr = (
        driver.Builder()
        .enable_dynamic_execution(allow_experimental_mode=True)
        .with_local_executor(executors.SynchronousLocalTaskExecutor())
        .with_remote_executor(executors.MultiThreadingExecutor(max_tasks=kwargs.get("global_parameters", {}).get("max_tasks", 1)))
        .with_modules(rbp_driver)
        .build()
    )
    dr.execute(
        ['close_rbp'],
        inputs={'config': kwargs, 'lookup_dict': lookup_dict, 'rbp_tables': rbp_tables}
    )


def run_tool(**kwargs):
    """Deprecated: Use run_annotation() instead."""
    run_annotation(**kwargs)

if __name__ == "__main__":
    cli()