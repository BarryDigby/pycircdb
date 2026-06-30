Inputs
=======

circRNA Formats
-----------------

pycircdb requires a 1-column text file with circRNAs in the following format:

.. code-block:: text

    chr1:1000-2000|+
    chr2:3000-4000
    chr3:5000-6000|-


Strand information is optional. Coordinate inputs are matched tolerantly: strand is honoured when supplied and ignored when absent, and each start coordinate is matched as-is as well as
+/-1 to absorb 0-based vs 1-based differences. You therefore do not need to declare the coordinate system of your input (0-based vs 1-based) - pycircdb will match your inputs to the database coordinates regardless of the coordinate system used.

Quantification Tool Cheatsheet
-----------------------------

It is worthwhile noting that different quantification tools report circRNAs as 0-based or 1-based coordinates.

pycircdb is 0-based. The following table may help reconcile your inputs vs. the outputs of pycircdb:

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
   * - circtools
     - 1-based
   * - DCC
     - 1-based
   * - find_circ
     - 0-based
   * - KNIFE
     - 1-based
   * - MapSplice
     - 0-based
   * - NCLscan
     - 1-based
   * - sailfish-cir
     - 1-based
   * - Segemehl
     - 0-based


.. caution::

    circRNAs on the + strand of "KNIFE","NCLscan" are reported as 'chr:end-start'.
    Please reverse the end and start coordinates to match the format of other tools.