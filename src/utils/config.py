"""
Configuration management utilities.
"""

import yaml
from pathlib import Path
from typing import Any, Dict, Optional
from functools import lru_cache


def load_config(config_path: str | Path) -> Dict[str, Any]:
    """
    Load a YAML configuration file.
    
    Parameters
    ----------
    config_path : str or Path
        Path to the YAML configuration file.
        
    Returns
    -------
    dict
        Configuration dictionary.
    """
    config_path = Path(config_path)
    
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    return config


@lru_cache(maxsize=8)
def get_config(config_name: str = "config") -> Dict[str, Any]:
    """
    Get a cached configuration by name.
    
    Parameters
    ----------
    config_name : str
        Name of the configuration file (without .yaml extension).
        Options: "config", "traits", "resources"
        
    Returns
    -------
    dict
        Configuration dictionary.
    """
    from src import CONFIG_DIR
    
    config_path = CONFIG_DIR / f"{config_name}.yaml"
    return load_config(config_path)


def get_trait_config(trait_name: str) -> Dict[str, Any]:
    """
    Get configuration for a specific trait.
    
    Parameters
    ----------
    trait_name : str
        Name of the trait (e.g., "LDL_cholesterol", "CAD").
        
    Returns
    -------
    dict
        Trait-specific configuration.
    """
    traits_config = get_config("traits")
    
    if trait_name not in traits_config.get("trait_tissue_priors", {}):
        raise ValueError(f"Unknown trait: {trait_name}")
    
    return traits_config["trait_tissue_priors"][trait_name]


def get_tissue_config(tissue_name: str) -> Dict[str, Any]:
    """
    Get configuration for a specific tissue.
    
    Parameters
    ----------
    tissue_name : str
        Name of the tissue (e.g., "liver", "adipose_visceral").
        
    Returns
    -------
    dict
        Tissue-specific configuration.
    """
    traits_config = get_config("traits")
    
    if tissue_name not in traits_config.get("tissue_ontology", {}):
        raise ValueError(f"Unknown tissue: {tissue_name}")
    
    return traits_config["tissue_ontology"][tissue_name]


def validate_config(config: Dict[str, Any]) -> bool:
    """
    Validate a configuration dictionary.
    
    Parameters
    ----------
    config : dict
        Configuration dictionary to validate.
        
    Returns
    -------
    bool
        True if valid, raises exception otherwise.
    """
    required_keys = ["project", "assembly", "directories", "gwas", "finemapping"]
    
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Missing required configuration key: {key}")
    
    # Validate assembly
    valid_assemblies = ["GRCh37", "GRCh38"]
    if config["assembly"] not in valid_assemblies:
        raise ValueError(f"Invalid assembly: {config['assembly']}. Must be one of {valid_assemblies}")
    
    return True
