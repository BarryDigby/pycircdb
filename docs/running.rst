Running pycircdb
================

Invokation
-----------

When running pycircdb on the command line, the order of the commands is important. The general syntax is:

.. code-block:: bash

    pycircdb [OPTIONS] COMMAND1 [ARGS]...


That is to say, global options (e.g. ``-c`` for the configuration file) must be specified before the command(s) (e.g. ``annotate``) and command-level options (e.g. ``-d`` for the database) must be specified directly after the command.

Multiple commands can be specified in a single invokation, for example:

.. code-block:: bash

    pycircdb -c config.json annotate -d 'circbase' fasta -d 'CSCD'


pycircdb Options
^^^^^^^^^^^^^^^^

View the set of options and commands available in pycircdb:

.. code-block:: bash

   pycircdb --help


Global Options
^^^^^^^^^^^^^^^^^^^

``-c, --config``
   Path to the JSON configuration file.

``-v, --verbose``
   Verbosity level. 0 = no output, 1 = errors only, 2 = errors and warnings,
   3 = all output (default=1).

``-h, --help``
   Show this help message and exit.

**Note that the configuration parameter is required by all commands except `init-demo`.**

Commands
^^^^^^^^^

``init-demo``
   Create a demo folder with example input data and a configuration file.

``annotate``
   Retrieve circRNA annotations from specified databases.

``fasta``
   Retrieve circRNA sequences in FASTA format from specified databases.

``mirna``
   Retrieve miRNA interactions for identified circRNAs.

``rbp``
   Retrieve RBP interactions for identified circRNAs.

Command-level Options
^^^^^^^^^^^^^^^^^^^^^

``init-demo``
    ``-f, --force``
        Overwrite the demo folder if it already exists.
    ``-h, --help``
        Show this help message and exit.
    
``annotate``
    ``-d, --database``
        Comma-separated list of databases to use for annotation. Options are
        ``arraystar``, ``circbank``, ``circbase``, ``circpedia``, ``circRNA_DB``, ``CSCD``, ``exorbase``
    ``-h, --help``
        Show this help message and exit.

``fasta``
    ``-d, --database``
        Comma-separated list of databases to use for sequence extraction. Options are
        ``arraystar``, ``circbank``, ``circbase``, ``circpedia``, ``circRNA_DB``, ``CSCD``
    ``-h, --help``
        Show this help message and exit.

``mirna``
    ``-a, --algorithm``
        Comma-separated list of algorithms to use for miRNA binding site prediction. Options are
        ``miRanda`` and/or ``PITA`` and/or ``TargetScan``
    ``-h, --help``
        Show this help message and exit.

``rbp``
    ``-h, --help``
        Show this help message and exit.
    


