Usage
=====

Once installed, you can run a minimal working demo of ``pycircdb`` using the following command:

.. code-block:: bash

   pycircdb init-demo
   pycircdb -c test_config.json -v 2 annotate -d 'arraystar,circbase' fasta -d 'arraystar,circbase' mirna -a 'miRanda,TargetScan' rbp


Resulting in the following output directory where the code was executed:

.. code-block:: bash

   results/
   └── demo
      └── demo_sample
         ├── arraystar.fasta
         ├── arraystar_hits.txt
         ├── circbase.fasta
         ├── circbase_hits.txt
         ├── hg38_chr1_mirna_hits.txt.gz
         └── hg38_chr1_rbp_hits.txt.gz

   2 directories, 6 files


Configuration
=============

To run ``pycircdb`` on your own data, you will need to create a configuration file (in JSON format) that specifies the paths to your input data and relevant global parameters.

.. code-block:: json

   {
      "global_parameters": {
         "max_tasks": 1,
         "output_dir": "results/",
         "tmp_dir": "tmp/"
      },
      "samples": {
         "sample_one": {
               "file_path": "path/to/sample1.txt",
               "reference": "hg19"
         },
         "sample_two": {
               "file_path": "path/to/sample2.txt",
               "reference": "hg38"
         }
      }
   }

        

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
               "reference": "hg19"
           }
       }
   }

The ``samples`` dictionary is required. Each sample points to a file of
circRNA coordinates and declares its genome ``reference`` (e.g. ``hg19``).
Coordinates are matched tolerantly, so the coordinate system (0- or 1-based)
and strand do not need to be declared.

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