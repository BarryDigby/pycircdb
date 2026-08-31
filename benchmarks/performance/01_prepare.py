#!/usr/bin/env python3
"""
01_prepare.py  —  Prepare inputs and configs for pycircdb performance benchmarking.

High-confidence pool (100 / 1,000 / 10,000):
  Source: raw BED files from glio_data/ (5-col: chrom start end strand count)
  Filter: per-tool read count >= MIN_READS AND detected by >= MIN_TOOLS tools.

Scaling pool (100k / 200k / 300k):
  Source: same raw BED files, chr filter only, no read-count or consensus filter.
  Sampled without replacement.
"""
import json
import random
import re
from collections import Counter
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

MIN_READS = 3   # per-tool read count threshold for high-confidence pool
MIN_TOOLS = 3   # tool consensus threshold for high-confidence pool
SEED      = 42

RAW_DIRS = {
    'SRR5133906': Path('test/SRR5133906'),
    'SRR5133907': Path('test/SRR5133907'),
}

OUT_DIR = Path('benchmarks/performance/inputs')
CFG_DIR = Path('benchmarks/performance/configs')
OUT_DIR.mkdir(parents=True, exist_ok=True)
CFG_DIR.mkdir(parents=True, exist_ok=True)

_CHR_RE = re.compile(r'^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$')

# ---------------------------------------------------------------------------
# High-confidence pool (>= MIN_READS per tool, >= MIN_TOOLS consensus)
# ---------------------------------------------------------------------------

hc_pool: set = set()

for sample, raw_dir in RAW_DIRS.items():
    if not raw_dir.exists():
        print(f'WARNING: {raw_dir} not found')
        continue
    tool_counts: Counter = Counter()
    for f in sorted(raw_dir.glob('*.bed')):
        for line in f.read_text().splitlines():
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            chrom, start, end, strand, count_str = parts[0], parts[1], parts[2], parts[3], parts[4]
            try:
                if int(count_str) < MIN_READS:
                    continue
            except ValueError:
                continue
            if not _CHR_RE.match(chrom):
                continue
            tool_counts[f'{chrom}:{start}-{end}|{strand}'] += 1
    for c, n in tool_counts.items():
        if n >= MIN_TOOLS:
            hc_pool.add(c)

# ---------------------------------------------------------------------------
# Scaling pool (all reads >= 1, chr filter only, from raw BED files)
# ---------------------------------------------------------------------------

scale_pool: set = set()

for sample, raw_dir in RAW_DIRS.items():
    if not raw_dir.exists():
        print(f"WARNING: {raw_dir} not found — scaling pool may be incomplete")
        continue
    for f in sorted(raw_dir.glob('*.bed')):
        for line in f.read_text().splitlines():
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            chrom, start, end, strand = parts[0], parts[1], parts[2], parts[3]
            if not _CHR_RE.match(chrom):
                continue
            scale_pool.add(f'{chrom}:{start}-{end}|{strand}')

hc_sorted    = sorted(hc_pool)
scale_sorted = sorted(scale_pool)

print(f'High-confidence pool (>={MIN_TOOLS} tools):  {len(hc_sorted):,}')
print(f'Scaling pool (all unique, chr-filtered):     {len(scale_sorted):,}')

# ---------------------------------------------------------------------------
# Create input files + configs
# ---------------------------------------------------------------------------

TIERS = [
    ('100',   100,      hc_sorted),
    ('1000',  1_000,    hc_sorted),
    ('10000', 10_000,   hc_sorted),
    ('100k',  100_000,  scale_sorted),
    ('200k',  200_000,  scale_sorted),
    ('300k',  300_000,  scale_sorted),
]

rng = random.Random(SEED)

for label, n, pool in TIERS:
    if n > len(pool):
        print(f'Skipping {label}: pool only has {len(pool):,}')
        continue
    sample = sorted(rng.sample(pool, n))
    input_path = OUT_DIR / f'bench_{label}.txt'
    input_path.write_text('\n'.join(sample) + '\n')
    cfg = {
        'global_parameters': {
            'max_tasks': 1,
            'output_dir': f'benchmarks/performance/results/bench_{label}',
            'tmp_dir': 'tmp',
        },
        'samples': {
            f'bench_{label}': {'file_path': str(input_path), 'reference': 'hg19'}
        }
    }
    (CFG_DIR / f'bench_{label}.json').write_text(json.dumps(cfg, indent=4))
    print(f'  bench_{label}.txt  ({n:,} circRNAs)')

print('\nDone. Run 02_benchmark.py next.')
