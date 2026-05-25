import json
import shutil
from pathlib import Path
import rich_click as click
from typing import Tuple, List, Optional
from config import create_config, load_config, DEFAULT_CONFIG_FILE

# Workflow stuff
from hamilton import driver
from hamilton.execution import executors

from utils.connect_s3.download_annotation_tables import fetch_annotation_tables
from utils.connect_s3.download_sequence_tables import fetch_sequence_tables
from utils.connect_s3.download_mirna_tables import fetch_mirna_tables

import utils.detect_inputs.detect_inputs_driver as instantiate_lookup_driver
import utils.annotate.annotate_driver as annotation_driver
import utils.fasta.sequence_driver as sequence_driver
import utils.mirna.mirna_driver as mirna_driver


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
def cli():
    """Main CLI tool."""
    pass

@cli.command('config')
@click.option(
    "-n",
    "--name",
    type=str,
    default="config",
    help="Name of the generated config file (without extension)."
)
def _init_create_config(name: str):
    """Interactive builder to initialize a configuration file."""
    # TODO: In the future, use `rich` or `questionary` here to prompt the user
    # for input files, reference choices (hg19/hg38), and zero_based logic.
    click.echo(f"Initializing interactive configuration builder for '{name}.json'...")
    create_config(name)


@cli.command('annotate')
@click.option(
    "-c",
    "--config",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    required=True,
    help="Path to the JSON config file containing workflow parameters."
)
@click.option(
    "--fasta",
    is_flag=True,
    help="Also output circRNA sequences in FASTA format."
)
@click.option(
    "--mirna",
    is_flag=True,
    help="Also output miRNA interactions in tabular format."
)
def annotate(config: str, fasta: bool, mirna: bool):
    """Annotate circRNAs using a JSON configuration file."""
    # 1. Load configuration purely from the provided JSON
    # TODO: Ensure types are correct. Validate the JSON
    cfg = load_config(config)

    # 2. Validate essential keys exist
    if not cfg.get('samples'):
        raise click.UsageError("Configuration file must contain a 'samples' dictionary.")

    # 3. Run annotation and get lookup results
    lookup_dict = run_annotation(**cfg)

    # 4. If --fasta flag, output FASTA using the same lookup results
    if fasta:
        run_fasta(lookup_dict=lookup_dict, **cfg)

    # 5. If --mirna flag, output miRNA interactions using the same lookup results
    if mirna:
        run_mirna(lookup_dict=lookup_dict, **cfg)


@cli.command('fasta')
@click.option(
    "-c",
    "--config",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    required=True,
    help="Path to the JSON config file containing workflow parameters."
)
@click.option(
    "--mirna",
    is_flag=True,
    help="Also output miRNA interactions in tabular format."
)
def fasta(config: str, mirna: bool):
    """Output circRNA sequences in FASTA format."""
    # 1. Load configuration purely from the provided JSON
    cfg = load_config(config)

    # 2. Validate essential keys exist
    if not cfg.get('samples'):
        raise click.UsageError("Configuration file must contain a 'samples' dictionary.")

    # 3. Run fasta output, capture lookup_dict for optional mirna step
    lookup_dict = run_fasta(**cfg)

    # 4. If --mirna flag, reuse the lookup results
    if mirna:
        run_mirna(lookup_dict=lookup_dict, **cfg)


@cli.command('mirna')
@click.option(
    "-c",
    "--config",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    required=True,
    help="Path to the JSON config file containing workflow parameters."
)
def mirna(config: str):
    """Output miRNA interactions for identified circRNAs."""
    cfg = load_config(config)

    if not cfg.get('samples'):
        raise click.UsageError("Configuration file must contain a 'samples' dictionary.")

    run_mirna(**cfg)

def run_annotation(**kwargs):
    """Run the annotation workflow."""
    print("Running annotation with configuration:")
    print(kwargs)
    for k, v in kwargs.items():
         print(f"  {k}: {v}")

    # Lookup tables
    lookup_dict = instantiate_lookup_driver.instantiate_driver(kwargs)

    # Pull annotation tables
    annotation_tables = fetch_annotation_tables(lookup_dict)

    # Annotate + write to file in parallel
    dr = (
        driver.Builder()
        .enable_dynamic_execution(allow_experimental_mode=True)
        .with_local_executor(executors.SynchronousLocalTaskExecutor())
        .with_remote_executor(executors.MultiThreadingExecutor(max_tasks=kwargs.get("cpus", 1)))
        .with_modules(annotation_driver)
        .build()
    )
    
    dr.execute(
        ['close_annotation'],
        inputs={'config': kwargs, 'annotation_tables': annotation_tables, 'lookup_results': lookup_dict}
    )
    
    return lookup_dict


def run_fasta(lookup_dict=None, **kwargs):
    """Generate FASTA output from circRNA sequences.
    
    Args:
        lookup_dict: Optional pre-computed lookup results. If None, will be generated from scratch.
        **kwargs: Configuration parameters.
    """
    print("Running FASTA generation with configuration:")
    print(kwargs)
    
    # Generate lookup tables if not provided (i.e., fasta running standalone)
    if lookup_dict is None:
        print("Computing lookup tables from scratch...")
        lookup_dict = instantiate_lookup_driver.instantiate_driver(kwargs)
    else:
        print("Using lookup tables from previous step...")
    
    sequence_tables = fetch_sequence_tables(lookup_dict)

    dr = (
        driver.Builder()
        .enable_dynamic_execution(allow_experimental_mode=True)
        .with_local_executor(executors.SynchronousLocalTaskExecutor())
        .with_remote_executor(executors.MultiThreadingExecutor(max_tasks=kwargs.get("cpus", 1)))
        .with_modules(sequence_driver)
        .build()
    )

    dr.execute(
        ['close_sequence'],
        inputs={'config': kwargs, 'lookup_dict': lookup_dict, 'sequence_tables': sequence_tables}
    )

    return lookup_dict


def run_mirna(lookup_dict=None, **kwargs):
    """Output miRNA interactions for identified circRNAs.

    Args:
        lookup_dict: Optional pre-computed lookup results. If None, will be generated from scratch.
        **kwargs: Configuration parameters.
    """
    print("Running miRNA interaction output with configuration:")
    print(kwargs)

    if lookup_dict is None:
        print("Computing lookup tables from scratch...")
        lookup_dict = instantiate_lookup_driver.instantiate_driver(kwargs)
    else:
        print("Using lookup tables from previous step...")

    mirna_tables = fetch_mirna_tables()

    dr = (
        driver.Builder()
        .enable_dynamic_execution(allow_experimental_mode=True)
        .with_local_executor(executors.SynchronousLocalTaskExecutor())
        .with_remote_executor(executors.MultiThreadingExecutor(max_tasks=kwargs.get("cpus", 1)))
        .with_modules(mirna_driver)
        .build()
    )
    dr.execute(
        ['close_mirna'],
        inputs={'config': kwargs, 'lookup_dict': lookup_dict, 'mirna_tables': mirna_tables}
    )


def run_tool(**kwargs):
    """Deprecated: Use run_annotation() instead."""
    run_annotation(**kwargs)

if __name__ == "__main__":
    cli()