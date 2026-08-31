#!/usr/bin/env python3
"""
02_benchmark.py  —  Run pycircdb on each benchmark tier and report results.

Steps:
  1. Pre-warm the S3 cache using the largest tier (300k) so all subsequent
     runs measure compute time only. -- SKIP_WARMUP=True to skip this step.
  2. Run each tier (100 / 1k / 10k / 100k / 200k / 300k) with all four
     commands (annotate fasta mirna rbp), capturing wall-clock time and
     peak RAM via /usr/bin/time -v.
  3. Count hits per tier:
       Annotation / FASTA : union of unique INPUT circRNAs with >= 1 hit
                            (fuzzy +-1bp matching back to input coordinates)
       miRNA / RBP        : unique hg38 circRNA coordinates in output files
  4. Write benchmarks/performance/results/summary.tsv.

Prerequisites:
  - Run 01_prepare.py first to generate inputs and configs.
  - tmp/ cache directory accessible (shared or local).
"""
import csv
import gzip
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, Set

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

TIERS = ['100', '1000', '10000', '100k', '200k', '300k']
CFG_DIR     = Path('benchmarks/performance/configs')
RESULTS_DIR = Path('benchmarks/performance/results')
INPUTS_DIR  = Path('benchmarks/performance/inputs')
SUMMARY     = RESULTS_DIR / 'summary.tsv'
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

PYCIRCDB_CMD = ['uv', 'run', 'main.py']
TIME_CMD     = '/usr/bin/time'
SKIP_WARMUP  = True   # set True to skip cache pre-warming

_COORD_RE = re.compile(r'^(.*?):(\d+)-(\d+)(\|[+-])?$')

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def fuzzy_keys(coord: str) -> set:
    m = _COORD_RE.match(coord)
    if not m:
        return set()
    ch, s, e, st = m.group(1), int(m.group(2)), m.group(3), m.group(4) or ''
    return {f'{ch}:{s+d}-{e}{st}' for d in (-1, 0, 1) if s + d >= 0}


def build_key_map(input_file: Path) -> Dict[str, Set[str]]:
    """fuzzy_key -> set of original input coords that generated it."""
    key_map: Dict[str, Set[str]] = {}
    for c in input_file.read_text().splitlines():
        c = c.strip()
        if c:
            for k in fuzzy_keys(c):
                key_map.setdefault(k, set()).add(c)
    return key_map


def _matched_set(key_map: Dict[str, Set[str]], hit_file: Path) -> Set[str]:
    if not hit_file.exists():
        return set()
    matched: Set[str] = set()
    with hit_file.open() as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            val = row.get('hg19', '')
            if val in key_map:
                matched.update(key_map[val])
            stripped = val.rsplit('|', 1)[0] if '|' in val else val
            if stripped in key_map:
                matched.update(key_map[stripped])
    return matched


def count_annotation_union(key_map: Dict[str, Set[str]], result_dir: Path) -> int:
    union: Set[str] = set()
    for f in sorted(result_dir.glob('*_hits.txt')):
        union.update(_matched_set(key_map, f))
    return len(union)


def count_fasta_union(key_map: Dict[str, Set[str]], result_dir: Path) -> int:
    FASTA_DBS = {
        'arraystar': 'arraystar',
        'circbank':  'circbank',
        'circbase':  'circbase',
        'circpedia': 'circpedia',
        'circRNADb': 'circRNADb',
        'cscd':      'hg38',     # CSCD FASTA headers are hg38 coords
    }
    union: Set[str] = set()
    for db, key_col in FASTA_DBS.items():
        hit_file = result_dir / f'{db}_hits.txt'
        for fasta_file in sorted(result_dir.glob(f'{db}_*.fasta')):
            seq_ids = {
                ln[1:].strip()
                for ln in fasta_file.read_text().splitlines()
                if ln.startswith('>')
            }
            if not seq_ids or not hit_file.exists():
                continue
            with hit_file.open() as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                for row in reader:
                    if row.get(key_col, '') in seq_ids:
                        val = row.get('hg19', '')
                        if val in key_map:
                            union.update(key_map[val])
                        stripped = val.rsplit('|', 1)[0] if '|' in val else val
                        if stripped in key_map:
                            union.update(key_map[stripped])
    return len(union)


def count_module_unique(result_dir: Path, glob: str) -> int:
    unique: Set[str] = set()
    for f in sorted(result_dir.glob(glob)):
        with gzip.open(f, 'rt') as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                c = row.get('circRNA', '')
                if c:
                    unique.add(c)
    return len(unique)


def parse_time_output(time_file: Path):
    """Return (wall_clock_s, peak_ram_mb) from /usr/bin/time -v output."""
    wall_s, ram_mb = None, None
    for line in time_file.read_text().splitlines():
        if 'Elapsed (wall clock)' in line:
            t = line.split()[-1]  # h:mm:ss or m:ss
            parts = t.split(':')
            if len(parts) == 2:
                wall_s = int(parts[0]) * 60 + float(parts[1])
            elif len(parts) == 3:
                wall_s = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
        if 'Maximum resident set size' in line:
            kb = int(line.split()[-1])
            ram_mb = round(kb / 1024, 1)
    return wall_s, ram_mb


def run_tier(label: str, warmup: bool = False) -> dict:
    cfg   = CFG_DIR / f'bench_{label}.json'
    infile_name = f'bench_{label}.txt'
    # map label to input filename
    infile = INPUTS_DIR / infile_name
    sample_dir = RESULTS_DIR / f'bench_{label}' / f'bench_{label}'
    time_file  = RESULTS_DIR / f'bench_{label}_time.txt'

    if not warmup:
        import shutil
        if (RESULTS_DIR / f'bench_{label}').exists():
            shutil.rmtree(RESULTS_DIR / f'bench_{label}')

    cmd = (
        [TIME_CMD, '-v', '-o', str(time_file)] +
        PYCIRCDB_CMD +
        ['-c', str(cfg), '-v', '0',
         'annotate', 'fasta', 'mirna', 'rbp']
    )
    print(f'  Running {"(warmup) " if warmup else ""}bench_{label} ...')
    result = subprocess.run(cmd, capture_output=True)
    if result.returncode not in (0, 1):
        print(f'  Warning: exit code {result.returncode}', file=sys.stderr)
    if result.stderr:
        # print last 20 lines of stderr to surface pycircdb warnings
        tail = result.stderr.decode(errors='replace').strip().splitlines()[-20:]
        for line in tail:
            print(f'  [stderr] {line}', file=sys.stderr)

    if warmup:
        return {}

    wall_s, ram_mb = parse_time_output(time_file)
    n_input = sum(1 for l in infile.read_text().splitlines() if l.strip())

    key_map = build_key_map(infile)
    ann  = count_annotation_union(key_map, sample_dir)
    fas  = count_fasta_union(key_map, sample_dir)
    mirna = count_module_unique(sample_dir, 'hg38_*_mirna_hits.txt.gz')
    rbp   = count_module_unique(sample_dir, 'hg38_*_rbp_hits.txt.gz')

    disk_kb = sum(
        f.stat().st_size for f in (RESULTS_DIR / f'bench_{label}').rglob('*') if f.is_file()
    ) // 1024
    disk_mb = round(disk_kb / 1024, 1)

    return {
        'Tier': label,
        'Input': n_input,
        'Annotation': ann,
        'FASTA': fas,
        'miRNA': mirna,
        'RBP': rbp,
        'Wall_clock_s': round(wall_s, 1) if wall_s else '',
        'Peak_RAM_GB': round(ram_mb / 1024, 2) if ram_mb else '',
        'Results_MB': disk_mb,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

print('=== Step 1: pre-warming S3 cache with 300k tier ===')
if SKIP_WARMUP:
    print('Skipped (SKIP_WARMUP=True).\n')
else:
    run_tier('300k', warmup=True)
    print('Cache warm.\n')

print('=== Step 2: benchmarking all tiers ===')
rows = []
for label in TIERS:
    row = run_tier(label)
    rows.append(row)
    print(
        f"  {label:<8} ann={row['Annotation']:,}  fas={row['FASTA']:,}  "
        f"mir={row['miRNA']:,}  rbp={row['RBP']:,}  "
        f"{row['Wall_clock_s']}s  {row['Peak_RAM_GB']}GB  {row['Results_MB']}MB"
    )

fields = ['Tier', 'Input', 'Annotation', 'FASTA', 'miRNA', 'RBP',
          'Wall_clock_s', 'Peak_RAM_GB', 'Results_MB']
with SUMMARY.open('w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader()
    w.writerows(rows)

print(f'\nSummary written to {SUMMARY}')
