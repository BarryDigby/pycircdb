# pycircdb

A command-line tool for circular RNA (circRNA) database analysis and network construction.

## Overview

`pycircdb` enables researchers to query and analyse circRNA interactions, building biological networks from user-supplied lists of circRNAs, miRNAs, genes, and RNA-binding proteins (RBPs).

## Features

- **circRNA annotation** — annotate circRNAs with genomic information
- **miRNA target analysis** — build circRNA–miRNA interaction networks
- **RNA-binding protein analysis** — build circRNA–RBP interaction networks
- **ceRNA network construction** — build competitive endogenous RNA (circRNA–miRNA–mRNA) networks

## Installation

```bash
pip install pycircdb
```

## Quick Start

```bash
# Annotate circRNAs
pycircdb --circrna circrnas.txt --annotate-circrnas --outdir results/

# circRNA–miRNA interaction analysis
pycircdb --circrna circrnas.txt --mirna mirnas.txt --mirna-targets --outdir results/

# RNA-binding protein analysis
pycircdb --circrna circrnas.txt --rbp rbps.txt --rna-binding-proteins --outdir results/

# Full ceRNA network
pycircdb --circrna circrnas.txt --mirna mirnas.txt --gene genes.txt --cerna-network --outdir results/

# Using a configuration file
pycircdb --config my_analysis.yaml
```

## Input Files

| Option | Description |
|---|---|
| `--circrna` | Path to a file containing circRNA identifiers (one per line) |
| `--mirna` | Path to a file containing miRNA identifiers (one per line) |
| `--gene` | Path to a file containing gene identifiers (one per line) |
| `--rbp` | Path to a file containing RBP identifiers (one per line) |

## Analysis Modules

| Option | Description | Required Input |
|---|---|---|
| `--annotate-circrnas` | Annotate circRNAs with genomic information | `--circrna` |
| `--mirna-targets` | Build circRNA–miRNA interaction network | `--circrna` or `--mirna` |
| `--rna-binding-proteins` | Build circRNA–RBP interaction network | `--circrna` or `--rbp` |
| `--cerna-network` | Build ceRNA network | At least one of `--circrna`, `--mirna`, `--gene` |

## Algorithm Filtering

### circRNA Detection Algorithms
```bash
# Filter circRNAs detected by specific algorithms
--circrna-algorithm ciri,circexplorer2
--circrna-set-logic AND   # circRNA must be detected by all algorithms (default)
--circrna-set-logic OR    # circRNA detected by any algorithm
```

Valid options: `circexplorer2`, `circrna_finder`, `find_circ`, `ciri`

### miRNA Prediction Algorithms
```bash
# Filter miRNA interactions by prediction algorithm
--mirna-algorithm miRanda,TargetScan
--mirna-set-logic AND           # default
--mirna-type 7mer-m8,8mer-1a   # TargetScan MRE site types
--mirna-mfe -25.0               # miRanda minimum free energy threshold (-62.0 to -0.41)
--mirna-score 150.0             # miRanda interaction score threshold (140.0 to 220.0)
```

Valid algorithms: `miRanda`, `TargetScan`
Valid MRE types: `6mer`, `7mer-m8`, `7mer-1a`, `8mer-1a`

### Gene/miRNA–mRNA Databases
```bash
--gene-database miRTarBase,TarBase
--gene-set-logic AND   # default
```

Valid databases: `DIANA`, `ElMMo`, `MicroCosm`, `PITA`, `PicTar`, `TarBase`, `TargetScan`, `miRDB`, `miRTarBase`, `miRanda`, `miRecords`

## Configuration File

All options can be specified in a YAML configuration file:

```yaml
# Input files
circrna_file: "data/circrnas.txt"
mirna_file: "data/mirnas.txt"
gene_file: "data/genes.txt"

# Analysis modules
annotate_circrnas: true
mirna_targets: true
cerna_network: true

# Algorithm filtering
circrna_algorithm: ["ciri", "circexplorer2"]
circrna_set_logic: "AND"
mirna_algorithm: ["miRanda", "TargetScan"]
mirna_mfe: -30.0

# Runtime settings
outdir: "results/"
workers: 8
verbose: 1
```

```bash
pycircdb --config my_analysis.yaml
```

See [CONFIG_GUIDE.md](CONFIG_GUIDE.md) for full configuration documentation.

## Runtime Options

| Option | Default | Description |
|---|---|---|
| `--outdir` | `output/` | Output directory for results |
| `--workers` | `2` | Number of parallel processes |
| `--verbose` | `0` | Verbosity level (0 = normal, 1 = verbose, 2 = debug) |
| `--quiet` | `false` | Suppress output except errors |
| `--config` | — | Path to YAML configuration file |
| `--save-config` | — | Save current settings to a YAML file |

## Python API

```python
from pycircdb import run
from pycircdb.config import PycircdbConfig

config = PycircdbConfig(
    circrna_file="data/circrnas.txt",
    mirna_file="data/mirnas.txt",
    mirna_targets=True,
    outdir="results/",
    workers=4,
)

run(config)
```

## License

See [LICENSE](LICENSE) for details.
