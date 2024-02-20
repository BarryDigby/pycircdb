#!/usr/bin/env python

""" config module. Holds a single copy of
config variables to be used across all other modules """
from typing import List, Dict, Optional, Union

import inspect

import logging
import os
import subprocess
import sys
from datetime import datetime

import importlib_metadata
import yaml
import pyaml_env
import pycircdb 

from pycircdb.utils.util_functions import strtobool
from importlib.metadata import version
from typing import List, Dict, Optional

version = version("pycircdb")

logger = logging.getLogger("pycircdb")

# Constants
PYCIRCDB_DIR = os.path.dirname(os.path.realpath(inspect.getfile(pycircdb)))

quiet: bool
output_dir = os.path.realpath(os.getcwd())


# Other variables that set only through the CLI
annotate = bool
data_dir: Optional[str] = output_dir
circrna_algorithm: List[str]
circrna_set_logic = str
workers: int


# Populating the variables above from the default pycircdb config
config_defaults_path = os.path.join(PYCIRCDB_DIR, "utils", "config_defaults.yaml")
with open(config_defaults_path) as f:
    _default_config = yaml.safe_load(f)
for c, v in _default_config.items():
    globals()[c] = v
