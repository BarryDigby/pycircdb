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
DEFAULT_CONFIG_FILE = CONFIG_DIR / "config.json"

class ToolConfig(TypedDict, total=False):
    input: List[str]
    reference: Union[str, List[str], None]
    zero_based: Union[bool, List[bool], None]
    cpu: int

def create_config(name: str):
    """Generate a default configuration file."""
    output_path = Path.cwd() / f"{name}.json"
    if DEFAULT_CONFIG_FILE.is_file():
        shutil.copy(DEFAULT_CONFIG_FILE, output_path)
        click.echo(f"Created configuration file at {output_path}")
    else:
        # Fallback if internal config.json is somehow missing
        default_data = {
            "global_parameters": {
                "cpu": 1,
                "output_dir": "results/"
            },
            "samples": {
                "sample_1": {
                    "file_path": "path/to/sample1.txt",
                    "reference": "hg38",
                    "zero_based": True
                }
            }
        }
        with open(output_path, "w") as f:
            json.dump(default_data, f, indent=4)
        click.echo(f"Created default configuration file at {output_path}")

def load_config(user_config_path: Optional[str] = None, verbose: int = 1) -> ToolConfig:
    """
    Load default config and override with a user config file if provided.
    """
    config: ToolConfig = {}

    # Load defaults
    if DEFAULT_CONFIG_FILE.is_file():
        with DEFAULT_CONFIG_FILE.open() as f:
            default_config = json.load(f) or {}
            config.update(default_config)

    # Load user config if provided
    if user_config_path:
        user_path = Path(user_config_path)
        if user_path.is_file():
            with user_path.open() as f:
                user_config = json.load(f) or {}
                config.update(user_config)

    # Verbosity check
    if verbose == 1:
        display_path = str(Path(user_config_path).resolve()) if user_config_path else 'defaults'
        console.print(Text(f"✓ Configuration File Loaded: {display_path}", style="bold green"))
    elif verbose >= 2:
        display_path = str(Path(user_config_path).resolve()) if user_config_path else 'defaults'
        console.print(Text(f"✓ Configuration File Loaded: {display_path}", style="bold green"))
        
        global_table = Table(show_header=True, header_style="bold cyan", box=box.ROUNDED, expand=True)
        global_params = config.get("global_parameters", {})
        for key in global_params.keys():
            global_table.add_column(str(key), justify="center", style="magenta")
        if global_params:
            global_table.add_row(*[str(val) for val in global_params.values()])
            
        sample_table = Table(show_header=True, header_style="bold blue", box=box.ROUNDED, expand=True)
        sample_table.add_column("Sample Name", style="bold green")
        sample_table.add_column("Input File", style="yellow")
        sample_table.add_column("Ref", style="cyan", justify="center")
        sample_table.add_column("0-based", style="magenta", justify="center")

        for sample_name, sample_info in config.get("samples", {}).items():
            sample_table.add_row(
                sample_name,
                str(sample_info.get("input", sample_info.get("file_path", ""))),
                str(sample_info.get("reference", "")),
                str(sample_info.get("zero_based", ""))
            )
            
        panel_group = Group(
            Text("Global Parameters:", style="bold white"),
            global_table,
            Text(""),
            Text("Samples:", style="bold white"),
            sample_table
        )

        console.print(Panel(panel_group, title="[bold white]Workflow Configuration[/bold white]", border_style="green", expand=False))

    return config