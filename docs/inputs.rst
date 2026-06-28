Inputs
=======

circRNA Formats
-----------------

pycircdb requires a 1-column text file with circRNAs in the following format:

.. code-block:: text

    chr1:1000-2000|+
    chr2:3000-4000|-
    chr3:5000-6000|+


Strand information is required, as all circRNA quantification tools report strand information.

Quantification Tool Cheatsheet
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Quantification Tool
     - Coordinate System
   * - CIRI2
     - 1-based
   * - CIRIquant
     - 1-based
   * - CircExplorer2
     - 0-based
   * - circRNA_finder
     - 0-based
   * - DCC
     - 1-based
   * - find_circ
     - 0-based
   * - MapSplice
     - 0-based
   * - Segemehl
     - 0-based

