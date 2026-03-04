from pathlib import Path
from typing import List, Optional, Union, Dict, Any
from dataclasses import dataclass, field
import yaml
import logging
from importlib.metadata import version

# Set up logging for the config module
logger = logging.getLogger(__name__)

class ConfigValidationError(Exception):
    """Custom exception for configuration validation errors"""
    pass

@dataclass
class PycircdbConfig:
    """
    Configuration class for pycircdb containing all settings and parameters.
    
    This class handles:
    - Default values for all configuration options
    - Loading/saving configuration from/to YAML files
    - Validation of parameter values and combinations
    - Type conversion and normalization
    """
    
    # Input files - paths to user data files
    circrna_file: Optional[Path] = None
    mirna_file: Optional[Path] = None  
    gene_file: Optional[Path] = None
    rbp_file: Optional[Path] = None
    
    # Module activation flags - which analyses to run
    annotate_circrnas: bool = False      # Annotate circRNAs with genomic information
    mirna_targets: bool = False          # Build circRNA-miRNA interaction networks
    rna_binding_proteins: bool = False   # Build circRNA-RBP interaction networks
    cerna_network: bool = False          # Build ceRNA (circRNA-miRNA-mRNA) networks
    
    # Filtering options for circRNA algorithms
    circrna_algorithm: Optional[List[str]] = None  # Which detection algorithms to include
    circrna_set_logic: str = "AND"                 # Logic for combining multiple algorithms
    
    # Filtering options for miRNA predictions
    mirna_algorithm: Optional[List[str]] = None    # Which prediction algorithms to use
    mirna_set_logic: str = "AND"                   # Logic for combining multiple algorithms
    mirna_type: Optional[List[str]] = None         # TargetScan MRE site types to include
    mirna_mfe: Optional[float] = None              # Minimum free energy threshold (miRanda)
    mirna_score: Optional[float] = None            # Interaction score threshold (miRanda)
    
    # Filtering options for gene databases
    gene_database: Optional[List[str]] = None      # Which miRNA-mRNA databases to use
    gene_set_logic: str = "AND"                    # Logic for combining multiple databases
    
    # Runtime configuration
    outdir: Path = Path("output")        # Output directory for results
    workers: int = 2                     # Number of parallel processes
    verbose: int = 0                     # Verbosity level (0=normal, 1=verbose, 2=debug)
    quiet: bool = False                  # Suppress most output except errors
    
    # Version information (set automatically)
    version: str = field(init=False)
    
    def __post_init__(self):
        """
        Called after __init__ to perform additional setup.
        Sets version and converts string paths to Path objects.
        """
        try:
            self.version = version("pycircdb")
        except Exception:
            self.version = "unknown"
        
        # Convert string paths to Path objects if needed
        for field_name in ['circrna_file', 'mirna_file', 'gene_file', 'rbp_file', 'outdir']:
            value = getattr(self, field_name)
            if value is not None and isinstance(value, str):
                setattr(self, field_name, Path(value))
    
    @classmethod
    def from_yaml(cls, config_file: Union[str, Path]) -> 'PycircdbConfig':
        """
        Load configuration from a YAML file.
        
        Args:
            config_file: Path to the YAML configuration file
            
        Returns:
            PycircdbConfig instance with loaded settings
            
        Raises:
            FileNotFoundError: If config file doesn't exist
            yaml.YAMLError: If YAML parsing fails
            ConfigValidationError: If loaded config is invalid
        """
        config_path = Path(config_file)
        
        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")
        
        logger.info(f"Loading configuration from: {config_path}")
        
        try:
            with open(config_path, 'r') as f:
                data = yaml.safe_load(f)
            
            # Handle empty or None YAML files before casting to dict
            if data is None:
                data = {}
            else:
                data = dict(data)  # Explicitly cast to dict for type checkers
            
            # Filter out None values and unknown keys
            valid_fields = {f.name for f in cls.__dataclass_fields__.values()}
            # Accept only known fields; set unknown fields to warning, allow None values
            filtered_data: Dict[str, Any] = {}
            for k, v in data.items():
                if k in valid_fields:
                    filtered_data[k] = v
                else:
                    logger.warning(f"Unknown config key in YAML: {k} (ignored)")
            
            # Create config instance
            config = cls(**filtered_data)
            
            # Validate the loaded configuration
            config.validate()
            
            return config
            
        except yaml.YAMLError as e:
            raise yaml.YAMLError(f"Error parsing YAML config file {config_path}: {e}")
    
    @classmethod
    def from_defaults(cls) -> 'PycircdbConfig':
        """
        Create configuration with default values loaded from config_defaults.yaml.
        
        Returns:
            PycircdbConfig instance with default settings
        """
        # Try to find the defaults file relative to this module
        defaults_path = Path(__file__).parent / "config_defaults.yaml"
        
        if defaults_path.exists():
            logger.debug(f"Loading defaults from: {defaults_path}")
            return cls.from_yaml(defaults_path)
        else:
            logger.warning("config_defaults.yaml not found, using hardcoded defaults")
            return cls()
    
    def to_yaml(self, config_file: Union[str, Path]) -> None:
        """
        Save current configuration to a YAML file.
        
        Args:
            config_file: Path where to save the configuration
            
        Raises:
            OSError: If file cannot be written
        """
        config_path = Path(config_file)
        
        # Create parent directories if they don't exist
        config_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Convert configuration to dictionary, handling special types
        data = {}
        for field_name, field_value in self.__dict__.items():
            if field_name == 'version':  # Skip version field in saved config
                continue
                
            if isinstance(field_value, Path):
                data[field_name] = str(field_value)
            elif field_value is None:
                data[field_name] = None
            else:
                data[field_name] = field_value
        
        logger.info(f"Saving configuration to: {config_path}")
        
        try:
            with open(config_path, 'w') as f:
                yaml.dump(data, f, default_flow_style=False, sort_keys=True)
        except OSError as e:
            raise OSError(f"Cannot write config file {config_path}: {e}")
    
    def validate(self) -> None:
        """
        Validate all configuration parameters.
        
        Raises:
            ConfigValidationError: If any validation check fails
        """
        errors: List[str] = []
        
        # Define valid options (these could be loaded from config_defaults.yaml)
        valid_circrna_algs = ['circexplorer2', 'circrna_finder', 'find_circ', 'ciri']
        valid_mirna_algs = ['miRanda', 'TargetScan']
        valid_mirna_types = ['6mer', '7mer-m8', '7mer-1a', '8mer-1a']
        valid_gene_dbs = [
            'DIANA', 'ElMMo', 'MicroCosm', 'PITA', 'PicTar', 'TarBase', 
            'TargetScan', 'miRDB', 'miRTarBase', 'miRanda', 'miRecords'
        ]
        valid_set_logic = ['AND', 'OR']
        
        # Validate input files exist if specified
        for file_field in ['circrna_file', 'mirna_file', 'gene_file', 'rbp_file']:
            file_path = getattr(self, file_field)
            if file_path is not None:
                if not isinstance(file_path, Path):
                    errors.append(f"{file_field} must be a Path object")
                elif not file_path.exists():
                    errors.append(f"{file_field} does not exist: {file_path}")
                elif not file_path.is_file():
                    errors.append(f"{file_field} is not a file: {file_path}")
        
        # Validate algorithm selections
        if self.circrna_algorithm is not None:
            for alg in self.circrna_algorithm:
                if alg not in valid_circrna_algs:
                    errors.append(f"Invalid circRNA algorithm: {alg}. Valid options: {valid_circrna_algs}")
        
        if self.mirna_algorithm is not None:
            for alg in self.mirna_algorithm:
                if alg not in valid_mirna_algs:
                    errors.append(f"Invalid miRNA algorithm: {alg}. Valid options: {valid_mirna_algs}")
        
        if self.mirna_type is not None:
            for mtype in self.mirna_type:
                if mtype not in valid_mirna_types:
                    errors.append(f"Invalid miRNA type: {mtype}. Valid options: {valid_mirna_types}")
        
        if self.gene_database is not None:
            for db in self.gene_database:
                if db not in valid_gene_dbs:
                    errors.append(f"Invalid gene database: {db}. Valid options: {valid_gene_dbs}")
        
        # Validate set logic options
        for logic_field in ['circrna_set_logic', 'mirna_set_logic', 'gene_set_logic']:
            logic_value = getattr(self, logic_field)
            if logic_value not in valid_set_logic:
                errors.append(f"Invalid {logic_field}: {logic_value}. Valid options: {valid_set_logic}")
        
        # Validate numeric ranges
        if self.mirna_mfe is not None:
            if not (-62.0 <= self.mirna_mfe <= -0.41):
                errors.append(f"mirna_mfe must be between -62.0 and -0.41, got: {self.mirna_mfe}")
        
        if self.mirna_score is not None:
            if not (140.0 <= self.mirna_score <= 220.0):
                errors.append(f"mirna_score must be between 140.0 and 220.0, got: {self.mirna_score}")
        
        # Validate workers is positive
        if self.workers <= 0:
            errors.append(f"workers must be positive, got: {self.workers}")
        
        # Validate verbose level
        if self.verbose < 0:
            errors.append(f"verbose level cannot be negative, got: {self.verbose}")
        
        # Validate module combinations and requirements
        self._validate_module_requirements(errors)
        
        # Raise exception if any errors found
        if errors:
            error_msg = "Configuration validation failed:\n" + "\n".join(f"  - {error}" for error in errors)
            raise ConfigValidationError(error_msg)
        
        logger.debug("Configuration validation passed")
    
    def _validate_module_requirements(self, errors: List[str]) -> None:
        """
        Validate that required input files are provided for selected modules.
        
        Args:
            errors: List to append validation errors to
        """
        # Check if any analysis modules are selected
        if not any([self.annotate_circrnas, self.mirna_targets, self.rna_binding_proteins, self.cerna_network]):
            errors.append("At least one analysis module must be selected")
        
        # Validate module-specific requirements
        if self.annotate_circrnas and self.circrna_file is None:
            errors.append("annotate_circrnas requires circrna_file to be specified")
        
        if self.mirna_targets:
            if self.circrna_file is None and self.mirna_file is None:
                errors.append("mirna_targets requires at least circrna_file or mirna_file")
        
        if self.rna_binding_proteins:
            if self.circrna_file is None and self.rbp_file is None:
                errors.append("rna_binding_proteins requires at least circrna_file or rbp_file")
        
        if self.cerna_network:
            if not any([self.circrna_file, self.mirna_file, self.gene_file]):
                errors.append("cerna_network requires at least one of: circrna_file, mirna_file, or gene_file")
    
    def update_from_cli(self, **kwargs: Any) -> None:
        """
        Update configuration with values from command line arguments.
        Only updates fields that are not None in kwargs.
        
        Args:
            **kwargs: Command line arguments to update
        """
        for key, value in kwargs.items():
            if value is not None and hasattr(self, key):
                # Convert string paths to Path objects
                if key.endswith('_file') or key == 'outdir':
                    if isinstance(value, str):
                        value = Path(value)
                
                setattr(self, key, value)
                logger.debug(f"Updated {key} from CLI: {value}")
    
    def get_summary(self) -> Dict[str, Any]:
        """
        Get a summary of current configuration for logging/reporting.
        
        Returns:
            Dictionary containing key configuration settings
        """
        return {
            'version': self.version,
            'input_files': {
                'circrna': str(self.circrna_file) if self.circrna_file else None,
                'mirna': str(self.mirna_file) if self.mirna_file else None,
                'gene': str(self.gene_file) if self.gene_file else None,
                'rbp': str(self.rbp_file) if self.rbp_file else None,
            },
            'modules': {
                'annotate_circrnas': self.annotate_circrnas,
                'mirna_targets': self.mirna_targets,
                'rna_binding_proteins': self.rna_binding_proteins,
                'cerna_network': self.cerna_network,
            },
            'filtering': {
                'circrna_algorithm': self.circrna_algorithm,
                'mirna_algorithm': self.mirna_algorithm,
                'gene_database': self.gene_database,
            },
            'runtime': {
                'outdir': str(self.outdir),
                'workers': self.workers,
                'verbose': self.verbose,
            }
        }