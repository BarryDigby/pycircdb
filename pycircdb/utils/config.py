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

from pycircdb.utils.util_functions import strtobool
from importlib.metadata import version
from typing import List, Dict, Optional

version = version("pycircdb")

logger = logging.getLogger("pycircdb")

quiet: bool
output_dir = os.path.realpath(os.getcwd())


# Other variables that set only through the CLI
data_dir: Optional[str] = output_dir