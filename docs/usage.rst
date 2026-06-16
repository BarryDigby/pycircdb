Usage
=====

``pycircdb`` is a command-line tool driven by a JSON configuration file. The
general invocation pattern is:

.. code-block:: bash

   pycircdb -c <config.json> <command> [options]

Subcommands can be chained together so that the input circRNAs are looked up
once and reused across multiple steps, for example:

.. code-block:: bash

   pycircdb -c config.json annotate fasta mirna rbp

Global options
--------------

``-c, --config``
   Path to the JSON config file containing workflow parameters. Required
   before any subcommand.

``-v, --verbose``
   Verbosity level: ``0`` (silent), ``1`` (high-level, default), ``2`` (all
   outputs).

``-h, --help``
   Show the help message and exit.

The configuration file
----------------------

The config file defines global parameters and the samples to process. A
minimal example:

.. code-block:: json

   {
       "global_parameters": {
           "max_tasks": 1,
           "output_dir": "results/ci_test/"
       },
       "samples": {
           "ci_test_sample": {
               "file_path": "test/tester.txt",
               "reference": "hg19",
               "zero_based": true
           }
       }
   }

The ``samples`` dictionary is required. Each sample points to a file of
circRNA coordinates and declares its genome ``reference`` (e.g. ``hg19``) and
whether the coordinates are ``zero_based``.

Commands
--------

annotate
~~~~~~~~

Annotate circRNAs using the configured databases.

.. code-block:: bash

   pycircdb -c config.json annotate -d arraystar,circbase,circatlas

``-d, --database``
   Comma-separated list of databases to use. Valid options: ``arraystar``,
   ``circatlas``, ``circbank``, ``circbase``, ``circpedia``, ``circrna_db``,
   ``cscd``, ``exorbase``.

fasta
~~~~~

Output circRNA sequences in FASTA format.

.. code-block:: bash

   pycircdb -c config.json fasta -d arraystar,circbase

``-d, --database``
   Comma-separated list of databases to use. Valid options: ``arraystar``,
   ``circatlas``, ``circbank``, ``circbase``, ``circpedia``, ``circrna_db``,
   ``cscd``.

mirna
~~~~~

Output miRNA interactions for the identified circRNAs.

.. code-block:: bash

   pycircdb -c config.json mirna -a miRanda,PITA,TargetScan

``-a, --algorithm``
   Comma-separated list of algorithms to use. Valid options: ``miRanda``,
   ``PITA``, ``TargetScan``.

rbp
~~~

Output RNA-binding protein (RBP) interactions for the identified circRNAs.

.. code-block:: bash

   pycircdb -c config.json rbp