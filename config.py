import json
import shutil
import rich_click as click
from datetime import datetime
from importlib.metadata import version as pkg_version
from rich.console import Console, Group
from rich.panel import Panel
from rich.text import Text
from rich.table import Table
from rich import box
from pathlib import Path
from typing import List, Optional, TypedDict, Union
from utils.md5sum_check import _get_db_prefix

db_version = _get_db_prefix()
db_version = db_version.replace("v", "").replace("_", ".")  # e.g. v1_0 -> 1.0

console = Console(stderr=True, highlight=False)
CONFIG_DIR = Path(__file__).parent.absolute()

DEFAULT_CONFIG_DATA = {
    "global_parameters": {
        "max_tasks": 1,
        "output_dir": "results/",
        "tmp_dir": "tmp"
    },
    "samples": {
        "sample_1": {
            "file_path": "path/to/sample1.txt",
            "reference": "hg38"
        }
    }
}

class ToolConfig(TypedDict, total=False):
    file_path: List[str]
    reference: Union[str, List[str], None]
    max_tasks: int

REQUIRED_SAMPLE_KEYS = ("file_path", "reference")

def _validate_samples(samples: dict) -> None:
    """
    Ensure every sample explicitly supplies the required keys.

    These are deliberately not defaulted: guessing the reference genome or
    coordinate system would silently produce zero lookup hits, so we fail loudly.
    """
    if not samples:
        raise click.UsageError("Configuration file must contain a 'samples' dictionary.")

    errors = []
    for sample_name, sample_info in samples.items():
        if not isinstance(sample_info, dict):
            errors.append(f"  - '{sample_name}': must be a mapping of sample parameters.")
            continue
        missing = [key for key in REQUIRED_SAMPLE_KEYS if sample_info.get(key) is None]
        if missing:
            errors.append(f"  - '{sample_name}': missing required key(s): {', '.join(missing)}.")

    if errors:
        raise click.UsageError(
            "Each sample must explicitly define "
            f"{', '.join(REQUIRED_SAMPLE_KEYS)}.\n" + "\n".join(errors)
        )

def create_config(name: str):
    """Generate a default configuration file."""
    output_path = Path.cwd() / f"{name}.json"
    with open(output_path, "w") as f:
        json.dump(DEFAULT_CONFIG_DATA, f, indent=4)
    click.echo(f"Created default configuration file at {output_path}")

def load_config(user_config_path: Optional[str] = None, verbose: int = 1) -> ToolConfig:
    """
    Load default config and override with a user config file if provided.
    """
    config: ToolConfig = {}

    # Load defaults directly from memory (copy so we never mutate the module-level defaults)
    config.update(json.loads(json.dumps(DEFAULT_CONFIG_DATA)))

    # Load user config if provided
    if user_config_path:
        user_path = Path(user_config_path)
        if user_path.is_file():
            with user_path.open() as f:
                user_config = json.load(f) or {}
                # Deep-merge global_parameters so a partial override keeps the remaining
                # defaults (e.g. supplying only max_tasks still yields default output_dir/tmp_dir).
                user_globals = user_config.pop("global_parameters", None)
                if user_globals is not None:
                    merged_globals = dict(config.get("global_parameters", {}))
                    merged_globals.update(user_globals)
                    config["global_parameters"] = merged_globals
                config.update(user_config)

    # Validate each sample, can't have missing keys here. 
    _validate_samples(config.get("samples", {}))

    # Verbosity check
    if verbose >= 1:
        display_path = str(Path(user_config_path).resolve()) if user_config_path else 'defaults'
        console.print(Text(f"✓ Configuration File Loaded: {display_path}", style="bold green"))

    return config

def print_config_panel(config: ToolConfig, user_config_path: Optional[str] = None):
    """Prints the rich panel for the workflow configuration."""
    global_table = Table(show_header=True, header_style="bold cyan", box=box.ROUNDED, expand=True)
    global_params = config.get("global_parameters", {})
    for key in global_params.keys():
        global_table.add_column(str(key), justify="center", style="magenta")
    if global_params:
        global_table.add_row(*[str(val) for val in global_params.values()])
        
    sample_table = Table(show_header=True, header_style="bold blue", box=box.ROUNDED, expand=True)
    sample_table.add_column("Sample Name", style="bold green")
    sample_table.add_column("File Path", style="yellow")
    sample_table.add_column("Ref", style="cyan", justify="center")

    for sample_name, sample_info in config.get("samples", {}).items():
        sample_table.add_row(
            sample_name,
            str(sample_info.get("file_path", "")),
            str(sample_info.get("reference", ""))
        )

    db_table = Table(show_header=True, header_style="bold yellow", box=box.ROUNDED, expand=True)
    db_table.add_column("Annotation DB", style="cyan")
    db_table.add_column("FASTA DB", style="magenta")
    db_table.add_column("miRNA Algorithms", style="green")

    ann_dbs = config.get("annotate_databases") or []
    fas_dbs = config.get("fasta_databases") or []
    mir_algs = config.get("mirna_algorithms") or []

    max_len = max(len(ann_dbs), len(fas_dbs), len(mir_algs))
    for i in range(max_len):
        a = ann_dbs[i] if i < len(ann_dbs) else ""
        f = fas_dbs[i] if i < len(fas_dbs) else ""
        m = mir_algs[i] if i < len(mir_algs) else ""
        db_table.add_row(a, f, m)
        
    panel_group = Group(
        Text("Global Parameters:", style="bold white"),
        global_table,
        Text(""),
        Text("Samples:", style="bold white"),
        sample_table,
        Text(""),
        Text("Databases & Algorithms:", style="bold white"),
        db_table,
        Text(f"Database version: {db_version}", style="bold white")
    )

    console.print(Panel(panel_group, title="[bold white]Workflow Configuration[/bold white]", border_style="green", expand=False))


# (genome_build, version, accessed)
_DB_META: dict[str, tuple[str, str, str]] = {
    "arraystar": ("hg19", "2.0", "January 2025"),
    "circatlas":  ("hg38", "3.0", "January 2026"),
    "circbank":   ("hg19", "1.0", "January 2026"),
    "circbase":   ("hg19", "1.0", "January 2026"),
    "circnet":    ("hg19", "2.0", "January 2026"),
    "circpedia":  ("hg38", "3.0", "October 2025"),
    "circrnadb":  ("hg19", "1.0", "March 2026"),
    "cscd":       ("hg38", "2.0", "January 2026"),
    "exorbase":   ("hg38", "1.0", "May 2026"),
}


def write_execution_report(config: ToolConfig, config_path: Optional[str] = None) -> None:
    """Write a plain-text execution report to the output directory."""
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    output_dir = Path(config.get("global_parameters", {}).get("output_dir", "results/"))
    if not output_dir.is_absolute():
        output_dir = Path.cwd() / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        tool_version = pkg_version("pycircdb")
    except Exception:
        tool_version = "unknown"

    lines = [
        "pycircdb execution report",
        "=" * 40,
        f"Timestamp        : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"pycircdb version : {tool_version}",
        f"Database version : {db_version}",
        f"Config file      : {Path(config_path).resolve() if config_path else 'defaults'}",
        "",
        "Global Parameters",
        "-" * 20,
    ]
    for k, v in config.get("global_parameters", {}).items():
        lines.append(f"  {k}: {v}")

    lines += ["", "Samples", "-" * 20]
    for name, info in config.get("samples", {}).items():
        lines.append(f"  {name}:")
        lines.append(f"    file_path : {info.get('file_path', '')}")
        lines.append(f"    reference : {info.get('reference', '')}")

    ann_dbs  = config.get("annotate_databases") or []
    fas_dbs  = config.get("fasta_databases")    or []
    mir_algs = config.get("mirna_algorithms")   or []
    run_rbp  = config.get("run_rbp", False)
    if ann_dbs or fas_dbs or mir_algs or run_rbp    :
        lines += ["", "Databases & Algorithms", "-" * 20]
        if ann_dbs:  lines.append(f"  Annotation databases : {', '.join(ann_dbs)}")
        if fas_dbs:  lines.append(f"  FASTA databases      : {', '.join(fas_dbs)}")
        if mir_algs: lines.append(f"  miRNA algorithms     : {', '.join(mir_algs)}")

    def _db_line(key: str) -> str:
        build, ver, accessed = _DB_META.get(key.lower(), ("unknown", "unknown", "unknown"))
        return f"    {key:<14}: {build}, version {ver}, accessed {accessed}"

    prov_sections: list[tuple[str, list[str]]] = []
    if ann_dbs:
        prov_sections.append(("Annotation", ann_dbs))
    if fas_dbs:
        prov_sections.append(("FASTA", fas_dbs))
    if mir_algs:
        prov_sections.append(("miRNA", ["circnet", "cscd"]))
    if run_rbp:
        prov_sections.append(("RBP", ["cscd"]))
    if prov_sections:
        lines += ["", "Database Provenance", "-" * 20]
        for section, dbs in prov_sections:
            lines.append(f"  {section}")
            for db in dbs:
                lines.append(_db_line(db))

    # ---------------------------------------------------------------------------
    # Mapping statistics — scanned from result files after all modules complete
    # ---------------------------------------------------------------------------
    import csv as _csv, gzip as _gz
    from collections import Counter as _Counter

    lines += ["", "Mapping Statistics", "-" * 30]

    for sample_name, sample_info in config.get("samples", {}).items():
        sample_dir = output_dir / sample_name
        if not sample_dir.exists():
            continue

        if len(config.get("samples", {})) > 1:
            lines.append(f"\n  [Sample: {sample_name}]")

        try:
            n_input = sum(1 for ln in Path(sample_info.get("file_path", "")).read_text().splitlines() if ln.strip())
        except Exception:
            n_input = None

        def _fmt(n: int) -> str:
            if n_input:
                return f"[{n}/{n_input}] {100 * n / n_input:.1f}%"
            return f"[{n} hits]"

        # Annotation: one hits file per database
        ann_files = sorted(sample_dir.glob("*_hits.txt"))
        if ann_files:
            lines.append("  Annotation")
            # circRNA -> native DB ID matching can be one-to-many in either build
            # (e.g. segmentally duplicated loci), which inflates raw hit counts
            # past 100%. Track duplicated hg19/hg38 loci per db so they can be
            # reported separately, and dedupe on the sample's own reference
            # build to keep the percentage meaningful.
            reference = str(sample_info.get("reference", "")).lower()
            multi_mappers: dict = {}
            for f in ann_files:
                db = f.stem.replace("_hits", "")
                with f.open(newline="") as fh:
                    rows = list(_csv.DictReader(fh, delimiter="\t"))
                fieldnames = rows[0].keys() if rows else []

                if reference in ("hg19", "hg38") and reference in fieldnames:
                    n = len({row[reference] for row in rows if row.get(reference, "").strip()})
                else:
                    n = len(rows)
                lines.append(f"    {db:<18} {_fmt(n)}")

                for build in ("hg19", "hg38"):
                    if build not in fieldnames:
                        continue
                    # Some source databases store a blank placeholder (e.g. a single
                    # whitespace char) where a coordinate failed to cross-map builds.
                    counts = _Counter(row[build] for row in rows if row.get(build, "").strip())
                    duplicates = sorted(val for val, cnt in counts.items() if cnt > 1)
                    if duplicates:
                        multi_mappers.setdefault(db, {})[build] = duplicates

            if multi_mappers:
                lines += ["", "    Multi-Mappers"]
                for db, builds in multi_mappers.items():
                    lines.append(f"      {db}")
                    for build, loci in builds.items():
                        lines.append(f"        {build}:")
                        for locus in loci:
                            lines.append(f"          - {locus}")

        # FASTA: one fasta file per database (e.g. arraystar_hg19.fasta)
        fasta_files = sorted(sample_dir.glob("*.fasta"))
        if fasta_files:
            lines += ["", "  FASTA"]
            for f in fasta_files:
                db = f.stem.rsplit("_", 1)[0]
                n = sum(1 for ln in f.open() if ln.startswith(">"))
                lines.append(f"    {db:<18} {_fmt(n)}")

        # miRNA: unique circRNAs per chromosome file
        mirna_files = sorted(sample_dir.glob("hg38_*_mirna_hits.txt.gz"))
        if mirna_files:
            lines += ["", "  miRNA"]
            for f in mirna_files:
                chrom = f.name.split("_mirna_")[0].replace("hg38_", "")
                try:
                    with _gz.open(f, "rt") as fh:
                        unique = {row.get("hg38", "") for row in _csv.DictReader(fh, delimiter="\t")}
                    n = len(unique - {""})
                except Exception:
                    n = 0
                lines.append(f"    {chrom:<18} {_fmt(n)}")

        # RBP: unique circRNAs per chromosome file
        rbp_files = sorted(sample_dir.glob("hg38_*_rbp_hits.txt.gz"))
        if rbp_files:
            lines += ["", "  RBP"]
            for f in rbp_files:
                chrom = f.name.split("_rbp_")[0].replace("hg38_", "")
                try:
                    with _gz.open(f, "rt") as fh:
                        unique = {row.get("hg38", "") for row in _csv.DictReader(fh, delimiter="\t")}
                    n = len(unique - {""})
                except Exception:
                    n = 0
                lines.append(f"    {chrom:<18} {_fmt(n)}")

    report_path = output_dir / f"execution_report_{timestamp}.txt"
    report_path.write_text("\n".join(lines) + "\n")
    console.print(Text(f"\u2713 Execution report written to {report_path}", style="bold green"))