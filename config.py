import json
import shutil
import rich_click as click
from rich.console import Console, Group
from rich.panel import Panel
from rich.text import Text
from rich.table import Table
from rich import box
from pathlib import Path
from typing import List, Optional, TypedDict, Union

console = Console(stderr=True, highlight=False)
CONFIG_DIR = Path(__file__).parent.absolute()

DEFAULT_CONFIG_DATA = {
    "global_parameters": {
        "max_tasks": 1,
        "output_dir": "results/",
        "tmp_dir": "tmp"
    },
    "samples": {
        "sample_1": {
            "file_path": "path/to/sample1.txt",
            "reference": "hg38",
            "zero_based": True
        }
    }
}

class ToolConfig(TypedDict, total=False):
    file_path: List[str]
    reference: Union[str, List[str], None]
    zero_based: Union[bool, List[bool], None]
    max_tasks: int

REQUIRED_SAMPLE_KEYS = ("file_path", "reference", "zero_based")

def _validate_samples(samples: dict) -> None:
    """
    Ensure every sample explicitly supplies the required keys.

    These are deliberately not defaulted: guessing the reference genome or
    coordinate system would silently produce zero lookup hits, so we fail loudly.
    """
    if not samples:
        raise click.UsageError("Configuration file must contain a 'samples' dictionary.")

    errors = []
    for sample_name, sample_info in samples.items():
        if not isinstance(sample_info, dict):
            errors.append(f"  - '{sample_name}': must be a mapping of sample parameters.")
            continue
        missing = [key for key in REQUIRED_SAMPLE_KEYS if sample_info.get(key) is None]
        if missing:
            errors.append(f"  - '{sample_name}': missing required key(s): {', '.join(missing)}.")

    if errors:
        raise click.UsageError(
            "Each sample must explicitly define "
            f"{', '.join(REQUIRED_SAMPLE_KEYS)}.\n" + "\n".join(errors)
        )

def create_config(name: str):
    """Generate a default configuration file."""
    output_path = Path.cwd() / f"{name}.json"
    with open(output_path, "w") as f:
        json.dump(DEFAULT_CONFIG_DATA, f, indent=4)
    click.echo(f"Created default configuration file at {output_path}")

def load_config(user_config_path: Optional[str] = None, verbose: int = 1) -> ToolConfig:
    """
    Load default config and override with a user config file if provided.
    """
    config: ToolConfig = {}

    # Load defaults directly from memory (copy so we never mutate the module-level defaults)
    config.update(json.loads(json.dumps(DEFAULT_CONFIG_DATA)))

    # Load user config if provided
    if user_config_path:
        user_path = Path(user_config_path)
        if user_path.is_file():
            with user_path.open() as f:
                user_config = json.load(f) or {}
                # Deep-merge global_parameters so a partial override keeps the remaining
                # defaults (e.g. supplying only max_tasks still yields default output_dir/tmp_dir).
                user_globals = user_config.pop("global_parameters", None)
                if user_globals is not None:
                    merged_globals = dict(config.get("global_parameters", {}))
                    merged_globals.update(user_globals)
                    config["global_parameters"] = merged_globals
                config.update(user_config)

    # Validate each sample, can't have missing keys here. 
    _validate_samples(config.get("samples", {}))

    # Verbosity check
    if verbose >= 1:
        display_path = str(Path(user_config_path).resolve()) if user_config_path else 'defaults'
        console.print(Text(f"✓ Configuration File Loaded: {display_path}", style="bold green"))

    return config

def print_config_panel(config: ToolConfig, user_config_path: Optional[str] = None):
    """Prints the rich panel for the workflow configuration."""
    global_table = Table(show_header=True, header_style="bold cyan", box=box.ROUNDED, expand=True)
    global_params = config.get("global_parameters", {})
    for key in global_params.keys():
        global_table.add_column(str(key), justify="center", style="magenta")
    if global_params:
        global_table.add_row(*[str(val) for val in global_params.values()])
        
    sample_table = Table(show_header=True, header_style="bold blue", box=box.ROUNDED, expand=True)
    sample_table.add_column("Sample Name", style="bold green")
    sample_table.add_column("File Path", style="yellow")
    sample_table.add_column("Ref", style="cyan", justify="center")
    sample_table.add_column("0-based", style="magenta", justify="center")

    for sample_name, sample_info in config.get("samples", {}).items():
        sample_table.add_row(
            sample_name,
            str(sample_info.get("file_path", "")),
            str(sample_info.get("reference", "")),
            str(sample_info.get("zero_based", ""))
        )

    db_table = Table(show_header=True, header_style="bold yellow", box=box.ROUNDED, expand=True)
    db_table.add_column("Annotation DB", style="cyan")
    db_table.add_column("FASTA DB", style="magenta")
    db_table.add_column("miRNA Algorithms", style="green")

    ann_dbs = config.get("annotate_databases", ["arraystar", "circbank", "circbase", "circpedia", "circrna_db", "cscd", "exorbase"])
    fas_dbs = config.get("fasta_databases", ["arraystar", "circbank", "circbase", "circpedia", "circrna_db", "cscd"])
    mir_algs = config.get("mirna_algorithms", ["miranda", "pita", "targetscan"])

    max_len = max(len(ann_dbs), len(fas_dbs), len(mir_algs))
    for i in range(max_len):
        a = ann_dbs[i] if i < len(ann_dbs) else ""
        f = fas_dbs[i] if i < len(fas_dbs) else ""
        m = mir_algs[i] if i < len(mir_algs) else ""
        db_table.add_row(a, f, m)
        
    panel_group = Group(
        Text("Global Parameters:", style="bold white"),
        global_table,
        Text(""),
        Text("Samples:", style="bold white"),
        sample_table,
        Text(""),
        Text("Databases & Algorithms:", style="bold white"),
        db_table
    )

    console.print(Panel(panel_group, title="[bold white]Workflow Configuration[/bold white]", border_style="green", expand=False))