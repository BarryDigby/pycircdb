Installation
============

You can install ``pycircdb`` using [uv](https://docs.astral.sh/uv/) (no separate Python installation required):

.. clode-block:: bash

   uv tool install pycircdb

Alternatively, you can install it from [PyPI](https://pypi.python.org/pypi/multiqc/) using ``pip`` (requires Python 3.10.1 or newer):
.. code-block:: bash

   pip install pycircdb


Requirements
=============

``pycircdb`` requires Python 3.10.1 or newer. It also requires the following Python packages:

.. code-block:: bash

   apache-hamilton>=1.90.0
   boto3>=1.42.89
   click>=8.3.2
   polars>=1.39.3
   rich-click>=1.9.7
