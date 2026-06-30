Configuration
==============

A JSON configuration file is required to run pycircdb using the CLI global option ``-c``.

At a high level, it specifies the output directory, the temporary directory and the input data for a workflow run. 

Global Parameters
^^^^^^^^^^^^^^^^^

``max_tasks``
    Maximum number of tasks to run in parallel. If not specified, the default is 1.
    
``output_dir``
   Path to the output directory. If it does not exist, it will be created. Defaults to ``./results``.

``tmp_dir``
   Path to the temporary directory. If it does not exist, it will be created. Defaults to ``./tmp``.

Note that the ``tmp_dir`` is designed to be re-used across multiple runs of pycircdb. 
This is because pycircdb will automatically scan the ``tmp_dir`` for any database files that have already been downloaded and will use these cached files for subsequent runs - facilitating offline runs in the future. 

Sample Parameters
^^^^^^^^^^^^^^^^^

``sample_name```
    Name of the sample being processed. This will be used to create a sub-directory in the output directory for the results of this run.

``input_path``
    Path to the input file containing circRNAs.

``reference``
    The reference genome used to perform circRNA quantification.
    

Please note that all of the sample parameters are required for each sample. If any of these parameters are missing, pycircdb will raise an error and exit.

Example Configuration File
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

    {
        "global_parameters": {
            "max_tasks": 4,
            "output_dir": "results/",
            "tmp_dir": "tmp/"
        },
        "samples": {
            "vromann": {
            "file_path": "test/vromann.txt",
            "reference": "hg38"
            },
            "glioblastoma_plus_dcc": {
            "file_path": "test/rnase_plus/glioblastoma_RNase_plus_dcc.txt",
            "reference": "hg38"
            },
            "glioblastoma_plus_ciriquant": {
            "file_path": "test/rnase_plus/glioblastoma_RNase_plus_ciriquant.txt",
            "reference": "hg38"
            }
        }
    }