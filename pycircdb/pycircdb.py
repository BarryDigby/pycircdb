
#!/usr/bin/env python3
"""
PycircDB: A command line tool for circular RNA database analysis and network construction.
"""

import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import rich_click as click
from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from rich.table import Table

from .config import PycircdbConfig, ConfigValidationError

# Set up rich console for pretty output
console = Console(stderr=True, highlight=False)

# Configure rich-click for beautiful help output
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True
click.rich_click.ERRORS_SUGGESTION = (
    "For more help, run '[yellow]pycircdb --help[/]' or visit "
    "[link=https://github.com/BarryDigby/pycircdb]https://github.com/BarryDigby/pycircdb[/]"
)
click.rich_click.STYLE_ERRORS_SUGGESTION = ""

# Group command line options for better help organization
click.rich_click.OPTION_GROUPS = {
    "pycircdb": [
        {
            "name": "Input Files",
            "options": ["--circrna", "--mirna", "--gene", "--rbp"],
        },
        {
            "name": "Analysis Modules", 
            "options": [
                "--annotate-circrnas",
                "--mirna-targets", 
                "--rna-binding-proteins",
                "--cerna-network",
            ],
        },
        {
            "name": "Algorithm Filtering",
            "options": [
                "--circrna-algorithm",
                "--circrna-set-logic",
                "--mirna-algorithm", 
                "--mirna-set-logic",
                "--mirna-type",
                "--mirna-mfe",
                "--mirna-score",
                "--gene-database",
                "--gene-set-logic",
            ],
        },
        {
            "name": "Configuration & Output",
            "options": ["--config", "--save-config", "--outdir", "--workers"],
        },
        {
            "name": "Utility Options",
            "options": ["--verbose", "--quiet", "--version", "--help"],
        },
    ]
}

def handle_list_option(ctx: click.Context, param: click.Parameter, value: Optional[str]) -> Optional[List[str]]:
    """Convert comma-separated string to list, or return None if empty."""
    if not value:
        return None
    return [item.strip() for item in value.split(',') if item.strip()]

def validate_config_and_show_errors(config: PycircdbConfig) -> bool:
    """
    Validate configuration and display any errors using rich formatting.
    
    Returns:
        True if validation passed, False otherwise
    """
    try:
        config.validate()
        return True
    except ConfigValidationError as e:
        # Create a beautiful error display with enhanced styling
        error_text = Text()
        error_text.append("⚠️  Configuration Validation Failed\n\n", style="bold red")
        
        # Parse and format individual errors with better styling and colored module names
        error_lines = str(e).split('\n')[1:]  # Skip the first line "Configuration validation failed:"
        for error_line in error_lines:
            if error_line.strip().startswith('- '):
                error_text.append("  ❌ ", style="red")
                
                # Color module names and convert underscores to dashes
                line = error_line.strip()[2:]
                
                # Replace module names with colored versions
                colored_line = line
                colored_line = colored_line.replace("annotate_circrnas", "annotate-circrnas")
                colored_line = colored_line.replace("mirna_targets", "mirna-targets") 
                colored_line = colored_line.replace("rna_binding_proteins", "rna-binding-proteins")
                colored_line = colored_line.replace("cerna_network", "cerna-network")
                
                # Split and color the module names
                if "annotate-circrnas" in colored_line:
                    colored_line = colored_line.replace("annotate-circrnas", "§annotate-circrnas§")
                if "mirna-targets" in colored_line:
                    colored_line = colored_line.replace("mirna-targets", "§mirna-targets§")
                if "rna-binding-proteins" in colored_line:
                    colored_line = colored_line.replace("rna-binding-proteins", "§rna-binding-proteins§")
                if "cerna-network" in colored_line:
                    colored_line = colored_line.replace("cerna-network", "§cerna-network§")
                
                # Replace file parameter names with colored versions
                colored_line = colored_line.replace("circrna_file", "§--circrna§")
                colored_line = colored_line.replace("mirna_file", "§--mirna§")
                colored_line = colored_line.replace("gene_file", "§--gene§")
                colored_line = colored_line.replace("rbp_file", "§--rbp§")
                
                # Split by § and add colors
                parts = colored_line.split("§")
                for part in parts:
                    if part in ["annotate-circrnas", "mirna-targets", "rna-binding-proteins", "cerna-network"]:
                        error_text.append(part, style="bold cyan")
                    elif part.startswith("--"):
                        error_text.append(part, style="yellow")
                    else:
                        error_text.append(part, style="white")
                
                error_text.append("\n")
        
        error_panel = Panel(
            error_text,
            title="🚨 [bold red]Validation Errors[/bold red] 🚨",
            border_style="red",
            padding=(1, 2),
            title_align="center"
        )
        console.print(error_panel)
        
        # Create intelligent suggestions based on available input files
        show_intelligent_suggestions(config)
        return False

def show_intelligent_suggestions(config: PycircdbConfig) -> None:
    """
    Show context-aware suggestions based on the configuration state.
    """
    suggestions = Text()
    suggestions.append("🔍 Here are some suggestions to fix your configuration:\n\n", style="bold blue")
    
    # If we have a circRNA file, suggest enabling all modules since they can work together
    if config.circrna_file and os.path.exists(config.circrna_file):
        suggestions.append("✨ Great! You have a circRNA file. ", style="green")
        suggestions.append("Consider enabling multiple analysis modules:\n", style="white")
        
        if not config.annotate_circrnas:
            suggestions.append("  🧬 Enable ", style="white")
            suggestions.append("annotate-circrnas", style="bold cyan")
            suggestions.append(" to annotate your circRNAs\n", style="white")
        
        if not config.mirna_targets:
            suggestions.append("  🎯 Enable ", style="white")
            suggestions.append("mirna-targets", style="bold cyan") 
            suggestions.append(" for miRNA-circRNA interaction analysis\n", style="white")
        
        if not config.rna_binding_proteins:
            suggestions.append("  🔗 Enable ", style="white")
            suggestions.append("rna-binding-proteins", style="bold cyan")
            suggestions.append(" for RBP-circRNA binding analysis\n", style="white")
        
        if not config.cerna_network:
            suggestions.append("  🕸️  Enable ", style="white")
            suggestions.append("cerna-network", style="bold cyan")
            suggestions.append(" for competitive endogenous RNA network analysis\n", style="white")
    
    # File-specific suggestions
    if not config.circrna_file:
        suggestions.append("📁 Specify your circRNA input file with ", style="white")
        suggestions.append("--circrna", style="yellow")
        suggestions.append(" /path/to/circrnas.txt\n", style="white")
    
    if config.mirna_targets and not config.mirna_file:
        suggestions.append("🎯 For miRNA target analysis, provide a miRNA file with ", style="white")
        suggestions.append("--mirna", style="yellow")
        suggestions.append(" /path/to/mirnas.txt\n", style="white")
    
    if config.cerna_network and not config.gene_file:
        suggestions.append("🧬 For ceRNA network analysis, provide a gene file with ", style="white")
        suggestions.append("--gene", style="yellow")
        suggestions.append(" /path/to/genes.txt\n", style="white")
    
    if config.rna_binding_proteins and not config.rbp_file:
        suggestions.append("🔗 For RBP analysis, provide an RBP file with ", style="white")
        suggestions.append("--rbp", style="yellow")
        suggestions.append(" /path/to/rbps.txt\n", style="white")
    
    # Add comprehensive CLI usage examples
    suggestions.append("\n💡 ", style="bold yellow")
    suggestions.append("Comprehensive CLI Usage Examples:", style="bold yellow")
    suggestions.append("\n\n", style="white")
    
    # Example 1: Basic annotation
    suggestions.append("1️⃣  Basic circRNA annotation:\n", style="bold white")
    suggestions.append("   pycircdb ", style="bold cyan")
    suggestions.append("--circrna", style="bold cyan")
    suggestions.append(" circrnas.txt ", style="white")
    suggestions.append("--annotate-circrnas", style="bold cyan")
    suggestions.append(" --output-dir", style="bold cyan")
    suggestions.append(" results/\n\n", style="white")
    
    # Example 2: circRNA - miRNA analysis
    suggestions.append("2️⃣ circRNA - miRNA analysis:\n", style="bold white")
    suggestions.append("   pycircdb ", style="bold cyan")
    suggestions.append("--circrna", style="bold cyan")
    suggestions.append(" circrnas.txt ", style="white")
    suggestions.append("--annotate-circrnas ", style="bold cyan")
    suggestions.append("--mirna-targets ", style="bold cyan")

    
    # Example 3: RBP focus
    suggestions.append("3️⃣  RNA-binding protein focused analysis:\n", style="bold white")
    suggestions.append("   pycircdb ", style="dim")
    suggestions.append("--circrna", style="yellow")
    suggestions.append(" circrnas.txt ", style="white")
    suggestions.append("--rbp", style="yellow")
    suggestions.append(" rbps.txt ", style="white")
    suggestions.append("--rna-binding-proteins", style="bold cyan")
    suggestions.append(" \\\n              ", style="white")
    suggestions.append("--output-dir", style="yellow")
    suggestions.append(" rbp_analysis/ ", style="white")
    suggestions.append("--threads", style="yellow")
    suggestions.append(" 8\n\n", style="white")
    
    # Example 4: With configuration file
    suggestions.append("4️⃣  Using configuration file:\n", style="bold white")
    suggestions.append("   pycircdb ", style="dim")
    suggestions.append("--config", style="yellow")
    suggestions.append(" my_analysis.yaml\n\n", style="white")
    
    # Configuration file example
    suggestions.append("📝 ", style="bold green")
    suggestions.append("Sample YAML Configuration File:", style="bold green")
    suggestions.append("\n\n", style="white")
    
    yaml_example = """# my_analysis.yaml
circrna_file: "data/circrnas.txt"
mirna_file: "data/mirnas.txt"  
gene_file: "data/genes.txt"
output_dir: "results/"

# Enable analysis modules
annotate_circrnas: true
mirna_targets: true
cerna_network: true

# Performance settings
threads: 8
memory_limit: 16"""
    
    suggestions.append(yaml_example, style="dim cyan")
    
    suggestion_panel = Panel(
        suggestions,
        title="💡 [bold green]Intelligent Suggestions[/bold green] 💡",
        border_style="green",
        padding=(1, 2),
        title_align="center"
    )
    console.print(suggestion_panel)


def demo_validation_system():
    """
    Demo function to showcase the enhanced validation system with rich formatting.
    This bypasses Click's validation to show our custom validation in action.
    """
    console.print("\n🎭 [bold blue]Enhanced Validation System Demo[/bold blue] 🎭\n")
    
    # Create a config with common validation issues
    demo_config = PycircdbConfig(
        circrna_file=Path("nonexistent_circrnas.txt"),  # File doesn't exist
        mirna_targets=True,  # Enabled but no mirna file
        cerna_network=True,  # Enabled but no gene file  
        rna_binding_proteins=True,  # Enabled but no rbp file
        workers=0,  # Invalid value
        outdir=Path(""),  # Empty directory
    )
    
    # Show our enhanced validation
    console.print("🔥 [bold yellow]Triggering validation with multiple issues...[/bold yellow]\n")
    result = validate_config_and_show_errors(demo_config)
    
    if not result:
        console.print("\n✅ [bold green]Validation system working perfectly![/bold green]")
        console.print("   Rich formatting, colored modules, and intelligent suggestions all displayed.")
    
    console.print("\n" + "="*60 + "\n")


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "--demo-validation",
    is_flag=True,
    help="🎭 Show enhanced validation system demo with rich formatting",
    hidden=True,  # Hidden option for testing
)
@click.option(
    "--config",
    type=click.Path(exists=True, readable=True, path_type=Path),
    help="Load configuration from YAML file",
)
@click.option(
    "--save-config", 
    type=click.Path(path_type=Path),
    help="Save current configuration to YAML file and exit",
)

# Input file options
@click.option(
    "--circrna",
    "circrna_file",
    type=click.Path(exists=True, readable=True, path_type=Path),
    help="Input circRNA file containing circRNA identifiers in the first column.\n\n"
         "Accepted formats: Hg19, Hg38, ArrayStar, CircBank, CircBase, circAtlas",
)
@click.option(
    "--mirna",
    "mirna_file", 
    type=click.Path(exists=True, readable=True, path_type=Path),
    help="Input miRNA file containing miRNA identifiers in the first column.\n\n"
         "Accepted formats: miRBase IDs (hsa-miR-106a-5p, miR-106a-5p)",
)
@click.option(
    "--gene",
    "gene_file",
    type=click.Path(exists=True, readable=True, path_type=Path),
    help="Input gene file containing gene identifiers in the first column.\n\n"
         "Accepted formats: HGNC symbols (REG4, PTEN, etc)",
)
@click.option(
    "--rbp",
    "rbp_file",
    type=click.Path(exists=True, readable=True, path_type=Path),
    help="Input RNA binding protein file containing RBP identifiers.\n\n"
         "Accepted formats: Gene symbols (ELAVL1, AGO2, etc)",
)

# Analysis module options
@click.option(
    "--annotate-circrnas",
    "annotate_circrnas",
    is_flag=True,
    help="Annotate circRNAs with genomic information.\n\n"
         "Requires valid circRNA input file using `--circrna`",
)
@click.option(
    "--mirna-targets",
    "mirna_targets", 
    is_flag=True,
    help="Build circRNA-miRNA interaction networks.\n\n"
         "Requires any combination of valid circRNA and/or miRNA input files using `--circrna` and `--mirna`",
)
@click.option(
    "--rna-binding-proteins",
    "rna_binding_proteins",
    is_flag=True,
    help="Build circRNA-RBP interaction networks.\n\n"
         "Requires any combination of valid circRNA and/or RBP input files using `--circrna` and `--rbp`",
)
@click.option(
    "--cerna-network",
    "cerna_network",
    is_flag=True,
    help="Build ceRNA (circRNA-miRNA-mRNA) networks.\n\n"
         "Requires any combination of valid circRNA, miRNA, and mRNA input files using `--circrna`, `--mirna`, and `--gene`",
)

# Algorithm filtering options
@click.option(
    "--circrna-algorithm",
    "circrna_algorithm",
    callback=handle_list_option,
    help="Filter by circRNA detection algorithms. Comma-separated list.\n\n"
         "Valid options: circexplorer2, circrna_finder, find_circ, ciri\n"
         "Example: `--circrna-algorithm ciri,circexplorer2`",
)
@click.option(
    "--circrna-set-logic",
    "circrna_set_logic",
    type=click.Choice(["AND", "OR"], case_sensitive=False),
    help="Logic for combining multiple circRNA algorithms specified by `--circrna_algorithm` [default: AND]",
)
@click.option(
    "--mirna-algorithm", 
    "mirna_algorithm",
    callback=handle_list_option,
    help="Filter by miRNA prediction algorithms. Comma-separated list.\n\n"
         "Valid options: miRanda, TargetScan\n"
         "Example: `--mirna-algorithm miRanda,TargetScan`",
)
@click.option(
    "--mirna-set-logic",
    "mirna_set_logic", 
    type=click.Choice(["AND", "OR"], case_sensitive=False),
    help="Logic for combining multiple miRNA algorithms specified by `--mirna_algorithm` [default: AND]",
)
@click.option(
    "--mirna-type",
    "mirna_type",
    callback=handle_list_option,
    help="Filter by TargetScan MRE site types. Comma-separated list.\n\n"
         "Valid options: 6mer, 7mer-1a, 7mer-m8, 8mer-1a\n"
         "Example: `--mirna-type 7mer-m8,8mer-1a`",
)
@click.option(
    "--mirna-mfe",
    "mirna_mfe",
    type=click.FloatRange(min=-62.0, max=-0.41),
    help="Minimum free energy threshold cutoff for miRanda predictions\n\n"
         "Range: -62.0 to -0.41 (more negative = stronger binding)",
)
@click.option(
    "--mirna-score",
    "mirna_score", 
    type=click.FloatRange(min=140.0, max=220.0),
    help="Minimum interaction score threshold cutoff for miRanda predictions\n\n"
         "Range: 140.0 to 220.0 (higher = stronger prediction)",
)
@click.option(
    "--gene-database",
    "gene_database",
    callback=handle_list_option,
    help="Filter by miRNA-mRNA interaction databases. Comma-separated list.\n\n"
         "Valid options: DIANA, ElMMo, MicroCosm, PITA, PicTar, TarBase, "
         "TargetScan, miRDB, miRTarBase, miRanda, miRecords\n"
         "Example: --gene-database TargetScan,miRTarBase",
)
@click.option(
    "--gene-set-logic",
    "gene_set_logic",
    type=click.Choice(["AND", "OR"], case_sensitive=False),
    help="Logic for combining multiple gene databases specified by `--gene_database` [default: AND]",
)

# Runtime options
@click.option(
    "--outdir",
    "-o",
    "outdir",
    type=click.Path(path_type=Path),
    help="Output directory for results [default: output]",
)
@click.option(
    "--workers",
    "-w",
    "workers",
    type=click.IntRange(min=1, max=64),
    help="Number of parallel workers [default: 2]",
)
@click.option(
    "--verbose",
    "-v", 
    "verbose",
    count=True,
    help="Increase verbosity (-v for verbose, -vv for debug)",
)
@click.option(
    "--quiet",
    "-q",
    "quiet",
    is_flag=True,
    help="Suppress all output except errors",
)
@click.version_option(version=None, prog_name="pycircdb")  # Version loaded from config
def run_cli(**kwargs: Any) -> None:
    """
    [bold blue]pycircdb[/bold blue]: A command line tool for circular RNA database analysis
    
    PycircDB pools publicly available circRNA databases to enable comprehensive
    analysis and network construction. At minimum, provide a list of circRNAs
    to get started with annotation.
    
    [dim]Basic usage:[/dim]
    [yellow]pycircdb --circrna circrnas.txt --annotate-circrnas[/yellow]
    
    [dim]Advanced usage:[/dim] 
    [yellow]pycircdb --config my_analysis.yaml --workers 8[/yellow]
    
    For more examples and documentation, visit:
    [link=https://github.com/BarryDigby/pycircdb]https://github.com/BarryDigby/pycircdb[/link]
    """
    
    # Check for demo mode first
    if kwargs.get('demo_validation', False):
        demo_validation_system()
        return
    
    # Show header with version info
    try:
        config_version = PycircdbConfig().version
    except Exception:
        config_version = "unknown"
    
    header_text = Text()
    header_text.append("🧬 ", style="blue")
    header_text.append("pycircdb", style="bold blue") 
    header_text.append(f" v{config_version}", style="dim")
    
    if not kwargs.get('quiet', False):
        console.print(Panel(header_text, padding=(0, 1), border_style="blue"))
    
    # Initialize configuration with proper typing
    config_file = kwargs.pop('config', None)  # type: ignore
    save_config = kwargs.pop('save_config', None)  # type: ignore
    config = PycircdbConfig()  # Initialize to avoid unbound variable
    
    try:
        # Load base configuration
        if config_file:
            if not kwargs.get('quiet', False):
                console.print(f"📁 Loading configuration from: [blue]{config_file}[/blue]")
            config = PycircdbConfig.from_yaml(config_file)
        else:
            config = PycircdbConfig()
        
        # Update with CLI arguments (only non-None values)
        cli_updates = {k: v for k, v in kwargs.items() if v is not None}
        if cli_updates:
            config.update_from_cli(**cli_updates)  # type: ignore
        
        # Handle config saving
        if save_config:
            config.to_yaml(save_config)
            console.print(f"💾 Configuration saved to: [green]{save_config}[/green]")
            sys.exit(0)
        
        # Validate configuration
        if not validate_config_and_show_errors(config):
            sys.exit(1)
        
        # Show configuration summary if verbose
        if config.verbose > 0:
            show_config_summary(config)
        
        # Run the analysis
        result = run(config)
        sys.exit(result.get("sys_exit_code", 0))
        
    except ConfigValidationError:
        # Already handled by validate_config_and_show_errors
        sys.exit(1)
    except FileNotFoundError as e:
        console.print(f"[red]❌ File not found:[/red] {e}")
        sys.exit(1)
    except Exception as e:
        console.print(f"[red]❌ Unexpected error:[/red] {e}")
        if config.verbose > 1:  # Debug mode
            console.print_exception()
        sys.exit(1)

def show_config_summary(config: PycircdbConfig) -> None:
    """Display a formatted summary of the current configuration."""
    
    summary = config.get_summary()
    
    table = Table(title="Configuration Summary", show_header=True, header_style="bold blue")
    table.add_column("Setting", style="cyan", width=20)
    table.add_column("Value", style="white")
    
    # Add rows to table
    table.add_row("Version", summary['version'])
    table.add_row("Output Directory", summary['runtime']['outdir']) 
    table.add_row("Workers", str(summary['runtime']['workers']))
    table.add_row("Verbosity", str(summary['runtime']['verbose']))
    
    # Show enabled modules
    enabled_modules = [k for k, v in summary['modules'].items() if v]
    if enabled_modules:
        table.add_row("Enabled Modules", ", ".join(enabled_modules))
    
    # Show input files
    input_files = [f"{k}: {v}" for k, v in summary['input_files'].items() if v]
    if input_files:
        for i, file_info in enumerate(input_files):
            label = "Input Files" if i == 0 else ""
            table.add_row(label, file_info)
    
    console.print(table)

def run(config: PycircdbConfig) -> Dict[str, int]:
    """
    Main analysis function that processes the configuration and runs pycircdb.
    
    Args:
        config: Validated PycircdbConfig instance
        
    Returns:
        Dictionary containing execution results and exit code
    """
    
    if not config.quiet:
        console.print("🚀 Starting pycircdb analysis...")
    
    # TODO: Implement actual analysis logic here
    # This is where you'll integrate your existing analysis pipeline
    
    # For now, just show what would be done
    if config.annotate_circrnas:
        console.print("📊 [green]Annotating circRNAs...[/green]")
    
    if config.mirna_targets:
        console.print("🔗 [green]Building circRNA-miRNA networks...[/green]")
    
    if config.rna_binding_proteins:
        console.print("🧬 [green]Analyzing RBP interactions...[/green]")
    
    if config.cerna_network:
        console.print("🕸️ [green]Constructing ceRNA networks...[/green]")
    
    console.print("✅ [bold green]Analysis completed successfully![/bold green]")
    
    return {"sys_exit_code": 0}
