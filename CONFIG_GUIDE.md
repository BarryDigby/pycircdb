---
output: pdf_document
---

# PycircDB Configuration System User Guide

## Overview

The enhanced PycircDB configuration system provides a robust, type-safe way to manage all application settings. It supports:

- **YAML configuration files** for reproducible analysis
- **Command-line overrides** for quick adjustments
- **Comprehensive validation** to catch errors early
- **Environment variable support** for deployment
- **Configuration saving/loading** for sharing settings

## Quick Start

### 1. Basic Usage

```python
from pycircdb.config import PycircdbConfig

# Create default configuration
config = PycircdbConfig()

# Create custom configuration
config = PycircdbConfig(
    workers=8,
    outdir=Path("my_results"),
    annotate_circrnas=True,
    circrna_algorithm=["ciri", "circexplorer2"]
)
```

### 2. Using Configuration Files

**Create a config file (`my_config.yaml`):**
```yaml
# Input files
circrna_file: "data/circrnas.txt"
mirna_file: "data/mirnas.txt"

# Analysis modules to run
annotate_circrnas: true
mirna_targets: true

# Algorithm filtering
circrna_algorithm: ["ciri", "circexplorer2"]
circrna_set_logic: "OR"
mirna_mfe: -25.0

# Runtime settings
outdir: "results_2024"
workers: 8
verbose: 1
```

**Load the configuration:**
```python
config = PycircdbConfig.from_yaml("my_config.yaml")
```

### 3. Command Line Integration

```python
import click
from pycircdb.config import PycircdbConfig

@click.command()
@click.option('--config', type=click.Path(exists=True), help='Config file')
@click.option('--workers', type=int, help='Number of workers')
@click.option('--outdir', type=str, help='Output directory')
# ... other options ...
def main(config=None, **kwargs):
    # Load base configuration
    if config:
        cfg = PycircdbConfig.from_yaml(config)
    else:
        cfg = PycircdbConfig()
    
    # Override with CLI arguments
    cfg.update_from_cli(**kwargs)
    
    # Validate before running
    cfg.validate()
    
    # Use configuration
    run_analysis(cfg)
```

## Configuration Reference

### Input Files
```yaml
circrna_file: null          # Path to circRNA identifiers file
mirna_file: null            # Path to miRNA identifiers file  
gene_file: null             # Path to gene identifiers file
rbp_file: null              # Path to RBP identifiers file
```

### Analysis Modules
```yaml
annotate_circrnas: false         # Annotate circRNAs with genomic info
mirna_targets: false             # Build circRNA-miRNA networks
rna_binding_proteins: false      # Build circRNA-RBP networks
cerna_network: false             # Build ceRNA networks
```

### Algorithm Filtering
```yaml
# CircRNA detection algorithms
circrna_algorithm: null          # ["circexplorer2", "circrna_finder", "find_circ", "ciri"]
circrna_set_logic: "AND"         # Logic: "AND" or "OR"

# MiRNA prediction algorithms  
mirna_algorithm: null            # ["miRanda", "TargetScan"]
mirna_set_logic: "AND"           # Logic: "AND" or "OR"
mirna_type: null                 # ["6mer", "7mer-1a", "7mer-m8", "8mer-1a"]
mirna_mfe: null                  # Range: -62.0 to -0.41 (miRanda only)
mirna_score: null                # Range: 140.0 to 220.0 (miRanda only)

# Gene/miRNA-mRNA databases
gene_database: null              # ["DIANA", "ElMMo", "MicroCosm", "PITA", "PicTar", 
                                 #  "TarBase", "TargetScan", "miRDB", "miRTarBase", 
                                 #  "miRanda", "miRecords"]
gene_set_logic: "AND"            # Logic: "AND" or "OR"
```

### Runtime Settings
```yaml
outdir: "output"            # Output directory path
workers: 2                  # Number of parallel processes
verbose: 0                  # Verbosity level (0=normal, 1=verbose, 2=debug)
quiet: false                # Suppress output except errors
```

## Validation Rules

The configuration system validates:

1. **File existence**: Input files must exist if specified
2. **Algorithm validity**: Only known algorithms are accepted
3. **Value ranges**: Numeric parameters must be within valid ranges
4. **Module requirements**: Required files must be provided for enabled modules
5. **Type safety**: All parameters must have correct types

### Module Requirements
- `annotate_circrnas`: requires `circrna_file`
- `mirna_targets`: requires `circrna_file` OR `mirna_file` 
- `rna_binding_proteins`: requires `circrna_file` OR `rbp_file`
- `cerna_network`: requires at least one of `circrna_file`, `mirna_file`, `gene_file`

## Advanced Usage

### Save Current Configuration
```python
config.to_yaml("my_saved_config.yaml")
```

### Load from Defaults
```python
config = PycircdbConfig.from_defaults()
```

### Get Configuration Summary
```python
summary = config.get_summary()
print(f"Version: {summary['version']}")
print(f"Enabled modules: {summary['modules']}")
```

### Handle Validation Errors
```python
from pycircdb.config import ConfigValidationError

try:
    config.validate()
except ConfigValidationError as e:
    print(f"Configuration errors:\n{e}")
    # Handle or fix errors
```

## Environment Variables

You can set configuration via environment variables using the `PYCIRCDB_` prefix:

```bash
export PYCIRCDB_WORKERS=16
export PYCIRCDB_VERBOSE=2
export PYCIRCDB_OUTDIR="/path/to/results"
```

## Migration from Old System

If you're upgrading from the old global configuration system:

**Old way:**
```python
import pycircdb.config as config
config.workers = 8
config.outdir = "results"
```

**New way:**
```python
from pycircdb.config import PycircdbConfig
config = PycircdbConfig(workers=8, outdir="results")
```

## Best Practices

1. **Use config files** for reproducible analyses
2. **Validate early** before starting long-running processes
3. **Save configurations** used for published results
4. **Use descriptive output directories** with timestamps/versions
5. **Version control your config files** alongside your data

## Example Workflows

### Interactive Analysis
```python
# Quick setup for exploration
config = PycircdbConfig(
    circrna_file="data.txt",
    annotate_circrnas=True,
    workers=4
)
```

### Production Pipeline
```yaml
# production_config.yaml
circrna_file: "/data/experiment_2024/circrnas.txt"
mirna_file: "/data/experiment_2024/mirnas.txt"
gene_file: "/data/experiment_2024/genes.txt"

annotate_circrnas: true
mirna_targets: true
cerna_network: true

circrna_algorithm: ["ciri", "circexplorer2"]
circrna_set_logic: "AND"
mirna_algorithm: ["miRanda", "TargetScan"]
mirna_mfe: -30.0

outdir: "/results/experiment_2024_v1"
workers: 32
verbose: 1
```

```bash
pycircdb --config production_config.yaml
```
