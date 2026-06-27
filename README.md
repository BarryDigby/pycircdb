<h1>
<picture>
  <source media="(prefers-color-scheme: light)" srcset="https://github.com/BarryDigby/pycircdb/raw/main/assets/pycircdb.png">
  <img src="https://github.com/BarryDigby/pycircdb/raw/main/assets/pycircdb.png" alt="pycircdb">
</picture>
</h1>

##### Explore the [documentation](https://pycircdb.readthedocs.io/en/latest/) to get up and running!

[![PyPI Version](https://img.shields.io/pypi/v/pycircdb)](https://pypi.python.org/pypi/pycircdb/)
[![DOI](https://img.shields.io/badge/DOI-pending-red.svg)](#)

---

Pycircdb is a python tool for annotating circular RNAs using four commands:

 - **annotation**: pycircdb queries Arraystar, Circbank, circBase, CIRCpedia, circRNADb, CSCD and exoRBase, generating outputs in the format returned by each database.
 - **fasta**: pycircdb queries Arraystar, CircBank, circBase, CIRCpedia, circRNADb and CSCD, generating outputs in FASTA format.
 - **mirna**: pycircdb queries CircNet and CSCD, generating bi-partite outputs of circRNAs-miRNAs.
 - **rbp**: pycircdb queries CSCD to return bi-partite outputs of circRNAs-RBPs. 

 ## Installation

You can install pycircdb using [uv](https://docs.astral.sh/uv/) (no separate Python installation required):

```bash
uv tool install pycircdb
```

Alternatively, install from [PyPI](https://pypi.python.org/pypi/multiqc/) using `pip`:

```bash
pip install pycircdb
```

## Usage

### Quickstart

Once installed, you can perform a minimal test-run to ensure the tool works as expected:

```bash
pycircdb init-demo
```

<picture>
  <source media="(prefers-color-scheme: light)" srcset="https://github.com/BarryDigby/pycircdb/raw/main/assets/pycircdb-init-demo.png">
  <img src="https://github.com/BarryDigby/pycircdb/raw/main/assets/pycircdb-init-demo.png" alt="pycircdb">
</picture>

Note that there are only three circRNAs in the demo input file. Whilst you can modify the parameters of the suggested command, it is possible that not all databases will return hits.  

```bash
pycircdb -c test_config.json -v 2 annotate -d 'arraystar,circbase' fasta -d 'arraystar,circbase' mirna -a 'miRanda,TargetScan' rbp
```

<picture>
  <source media="(prefers-color-scheme: light)" srcset="https://github.com/BarryDigby/pycircdb/raw/main/assets/pycircdb-demo-run.png">
  <img src="https://github.com/BarryDigby/pycircdb/raw/main/assets/pycircdb-demo-run.png" alt="pycircdb">
</picture>


The demo will produce the following output directory:

```bash
/home/barry/my_analysis/results/
└── demo
    └── demo_sample
        ├── arraystar.fasta
        ├── arraystar_hits.txt
        ├── circbase.fasta
        ├── circbase_hits.txt
        ├── hg38_chr1_mirna_hits.txt.gz
        └── hg38_chr1_rbp_hits.txt.gz

2 directories, 6 files
```

### Configuration

pycircdb requires an input configuration JSON file in order to locate input samples and define global and per-sample configuration settings.

Please see the [documentation](https://pycircdb.readthedocs.io/en/latest/) for a full description of the config file.

```json
{
  "global_parameters": {
    "max_tasks": 4,
    "output_dir": "results/",
    "tmp_dir": "tmp/"
  },
  "samples": {
    "vromann": {
      "file_path": "test/vromann.txt",
      "reference": "hg38",
      "zero_based": true
    },
    "glioblastoma_plus_dcc": {
      "file_path": "test/rnase_plus/glioblastoma_RNase_plus_dcc.txt",
      "reference": "hg38",
      "zero_based": true
    },
    "glioblastoma_plus_ciriquant": {
      "file_path": "test/rnase_plus/glioblastoma_RNase_plus_ciriquant.txt",
      "reference": "hg38",
      "zero_based": true
    }
  }
}
```

### Commands

For a full description of pycircdb commands and their options, please refer to the [documentation](https://pycircdb.readthedocs.io/en/latest/).


## Citation

Please consider citing pycircdb if you use it in your analysis

> **pycircdb: integrated circRNA database annotation for computational workflows** <br> _Barry Digby, Stephen Finn, Pilib Ó Broin_ <br>
> Pending (2026) <br>
> doi: [pending](pending)) <br>
> PMID: [pending](pending)

```BibTeX
@article{pending,
 author = {Digby, Barry and Finn, Stephen and Ó Broin, Pilib},
 title = {pycircdb: integrated circRNA database annotation for computational workflows},
 journal = {pending},
 volume = {pending},
 number = {pending},
 pages = {pending},
 year = {pending},
 doi = {pending},
 URL = {pending},
 eprint = {pending}
}
```

## Contributions & Support

Contributions, particularly new circRNA database resources are welcome, as are bug reports!
Please create a new [issue](https://github.com/BarryDigby/pycircdb/issues) for any of these, preferably with the `-v 2` flag enabled for richer logs.
Pull-requests for fixes and additions are very welcome.

This work was funded by Science Foundation Ireland, grant number 18/CRT/6214.