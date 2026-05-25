import json
import shutil
import rich_click as click
from pathlib import Path
from typing import List, Optional, TypedDict, Union

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

def load_config(user_config_path: Optional[str] = None) -> ToolConfig:
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

    return config