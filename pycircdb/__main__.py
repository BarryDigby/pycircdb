#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" __main__.py
~~~~~~~~~~~~~~~~~~~~
Called when pycircdb namespace is used, primarily by the cli:

$ pycricdb .
$ python -m pycircdb .
"""


from importlib_metadata import entry_points

from . import pycircdb


def run_pycircdb():
    # Call the main function
    pycircdb.run_cli(prog_name="pycircdb")


# Script is run directly
if __name__ == "__main__":
    run_pycircdb()
