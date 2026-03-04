#!/usr/bin/env python3
"""
Example script demonstrating how to use the enhanced PycircdbConfig class.

This shows various ways to create, validate, load, and save configurations.
"""

from pathlib import Path
from pycircdb.config import PycircdbConfig, ConfigValidationError

def main():
    """Demonstrate various configuration usage patterns."""
    
    print("=== PycircdbConfig Usage Examples ===\n")
    
    # Initialize variables to avoid "possibly unbound" warnings
    custom_config = None
    config_file = Path("example_config.yaml")
    
    # Example 1: Create default configuration
    print("1. Creating default configuration:")
    try:
        config = PycircdbConfig()
        print(f"   Version: {config.version}")
        print(f"   Default output dir: {config.outdir}")
        print(f"   Default workers: {config.workers}")
        print("   ✓ Default config created successfully\n")
    except Exception as e:
        print(f"   ✗ Error: {e}\n")
    
    # Example 2: Create configuration with custom values
    print("2. Creating custom configuration:")
    try:
        custom_config = PycircdbConfig(
            outdir=Path("custom_results"),
            workers=8,
            annotate_circrnas=True,
            mirna_targets=True,
            circrna_algorithm=["ciri", "circexplorer2"],
            circrna_set_logic="OR",
            verbose=1
        )
        print(f"   Output dir: {custom_config.outdir}")
        print(f"   Workers: {custom_config.workers}")
        print(f"   Modules enabled: annotate_circrnas={custom_config.annotate_circrnas}, mirna_targets={custom_config.mirna_targets}")
        print("   ✓ Custom config created successfully\n")
    except Exception as e:
        print(f"   ✗ Error: {e}\n")
    
    # Example 3: Save configuration to YAML
    print("3. Saving configuration to YAML:")
    try:
        if custom_config is not None:
            custom_config.to_yaml(config_file)
            print(f"   ✓ Configuration saved to: {config_file}")
            
            # Show the saved content
            with open(config_file) as f:
                content = f.read()
            print("   Content preview:")
            for line in content.split('\n')[:10]:  # Show first 10 lines
                print(f"     {line}")
            print("     ...\n")
        else:
            print("   ✗ No config to save\n")
    except Exception as e:
        print(f"   ✗ Error: {e}\n")
    
    # Example 4: Load configuration from YAML
    print("4. Loading configuration from YAML:")
    try:
        if config_file.exists():
            loaded_config = PycircdbConfig.from_yaml(config_file)
            print(f"   Loaded output dir: {loaded_config.outdir}")
            print(f"   Loaded workers: {loaded_config.workers}")
            print("   ✓ Configuration loaded successfully\n")
        else:
            print("   ✗ Config file not found\n")
    except Exception as e:
        print(f"   ✗ Error: {e}\n")
    
    # Example 5: Load from defaults file
    print("5. Loading from defaults file:")
    try:
        defaults_config = PycircdbConfig.from_defaults()
        print(f"   Default output dir: {defaults_config.outdir}")
        print("   ✓ Defaults loaded successfully\n")
    except Exception as e:
        print(f"   ✗ Error: {e}\n")
    
    # Example 6: Validation - valid config
    print("6. Validating valid configuration:")
    try:
        valid_config = PycircdbConfig(
            annotate_circrnas=True,
            circrna_algorithm=["ciri"],
            workers=4,
            mirna_mfe=-25.0,
            mirna_score=150.0
        )
        valid_config.validate()
        print("   ✓ Configuration validation passed\n")
    except ConfigValidationError as e:
        print(f"   ✗ Validation error: {e}\n")
    except Exception as e:
        print(f"   ✗ Unexpected error: {e}\n")
    
    # Example 7: Validation - invalid config
    print("7. Validating invalid configuration:")
    try:
        invalid_config = PycircdbConfig(
            annotate_circrnas=True,
            circrna_algorithm=["invalid_algorithm"],  # Invalid algorithm
            workers=-1,  # Invalid worker count
            mirna_mfe=-100.0,  # Out of range
            mirna_score=50.0   # Out of range
        )
        invalid_config.validate()
        print("   ✗ Validation should have failed!\n")
    except ConfigValidationError as e:
        print("   ✓ Validation correctly caught errors:")
        for line in str(e).split('\n')[1:6]:  # Show first few errors
            print(f"     {line}")
        print("   ...\n")
    except Exception as e:
        print(f"   ✗ Unexpected error: {e}\n")
    
    # Example 8: Update from CLI arguments
    print("8. Updating config from CLI arguments:")
    try:
        base_config = PycircdbConfig()
        print(f"   Before: workers={base_config.workers}, verbose={base_config.verbose}")
        
        # Simulate CLI arguments with explicit types
        base_config.update_from_cli(
            workers=16,
            verbose=2,
            outdir='cli_output',
            annotate_circrnas=True
        )
        print(f"   After: workers={base_config.workers}, verbose={base_config.verbose}")
        print(f"   Output dir: {base_config.outdir}")
        print("   ✓ CLI update successful\n")
    except Exception as e:
        print(f"   ✗ Error: {e}\n")
    
    # Example 9: Get configuration summary
    print("9. Getting configuration summary:")
    try:
        if custom_config is not None:
            summary = custom_config.get_summary()
            print("   Configuration summary:")
            print(f"     Version: {summary['version']}")
            print(f"     Modules: {summary['modules']}")
            print(f"     Runtime: {summary['runtime']}")
            print("   ✓ Summary generated successfully\n")
        else:
            print("   ✗ No config available for summary\n")
    except Exception as e:
        print(f"   ✗ Error: {e}\n")
    
    # Clean up
    try:
        if config_file.exists():
            config_file.unlink()
            print("✓ Cleaned up example files")
    except Exception:
        pass

if __name__ == "__main__":
    main()
