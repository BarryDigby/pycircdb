import logging
import os
from importlib.metadata import version
from typing import List, Dict, Optional

version = version("pycircdb")

logger = logging.getLogger("pycircdb")

quiet: bool
output_dir = os.path.realpath(os.getcwd())


# Other variables that set only through the CLI
data_dir: Optional[str] = None