"""
Input/Output utilities for genomic data.
"""

import gzip
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd
import numpy as np


def read_sumstats(
    filepath: str | Path,
    sep: str = "\t",
    compression: str = "infer",
    usecols: Optional[List[str]] = None,
    dtype: Optional[Dict[str, Any]] = None,
    nrows: Optional[int] = None,
) -> pd.DataFrame:
    """
    Read GWAS summary statistics file.
    
    Parameters
    ----------
    filepath : str or Path
        Path to the summary statistics file.
    sep : str
        Column separator.
    compression : str
        Compression type ('gzip', 'infer', None).
    usecols : list, optional
        Columns to read.
    dtype : dict, optional
        Column data types.
    nrows : int, optional
        Number of rows to read (for testing).
        
    Returns
    -------
    pd.DataFrame
        Summary statistics dataframe.
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    
    # Default dtypes for common columns
    default_dtype = {
        "chr": str,
        "chromosome": str,
        "CHR": str,
        "pos": "Int64",
        "position": "Int64",
        "BP": "Int64",
        "rsid": str,
        "SNP": str,
        "variant_id": str,
    }
    
    if dtype:
        default_dtype.update(dtype)
    
    df = pd.read_csv(
        filepath,
        sep=sep,
        compression=compression,
        usecols=usecols,
        dtype=default_dtype,
        nrows=nrows,
        low_memory=False,
    )
    
    return df


def write_sumstats(
    df: pd.DataFrame,
    filepath: str | Path,
    sep: str = "\t",
    compress: bool = True,
    index: bool = False,
) -> None:
    """
    Write GWAS summary statistics to file.
    
    Parameters
    ----------
    df : pd.DataFrame
        Summary statistics dataframe.
    filepath : str or Path
        Output file path.
    sep : str
        Column separator.
    compress : bool
        Whether to gzip compress the output.
    index : bool
        Whether to write row index.
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    
    if compress and not str(filepath).endswith('.gz'):
        filepath = Path(str(filepath) + '.gz')
    
    compression = 'gzip' if compress else None
    
    df.to_csv(
        filepath,
        sep=sep,
        compression=compression,
        index=index,
    )


def read_bed(
    filepath: str | Path,
    names: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Read BED format file.
    
    Parameters
    ----------
    filepath : str or Path
        Path to BED file.
    names : list, optional
        Column names.
        
    Returns
    -------
    pd.DataFrame
        BED data as dataframe.
    """
    filepath = Path(filepath)
    
    if names is None:
        names = ["chrom", "start", "end", "name", "score", "strand"]
    
    # Detect compression
    compression = "gzip" if str(filepath).endswith(".gz") else None
    
    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        names=names[:6] if len(names) > 6 else names,
        compression=compression,
        comment="#",
    )
    
    return df


def write_json(
    data: Union[Dict, List],
    filepath: str | Path,
    indent: int = 2,
    compress: bool = False,
) -> None:
    """
    Write data to JSON file.
    
    Parameters
    ----------
    data : dict or list
        Data to write.
    filepath : str or Path
        Output file path.
    indent : int
        JSON indentation level.
    compress : bool
        Whether to gzip compress.
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    
    # Custom JSON encoder for numpy types
    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)
    
    json_str = json.dumps(data, indent=indent, cls=NumpyEncoder)
    
    if compress:
        if not str(filepath).endswith('.gz'):
            filepath = Path(str(filepath) + '.gz')
        with gzip.open(filepath, 'wt') as f:
            f.write(json_str)
    else:
        with open(filepath, 'w') as f:
            f.write(json_str)


def read_json(filepath: str | Path) -> Union[Dict, List]:
    """
    Read JSON file.
    
    Parameters
    ----------
    filepath : str or Path
        Path to JSON file.
        
    Returns
    -------
    dict or list
        Parsed JSON data.
    """
    filepath = Path(filepath)
    
    if str(filepath).endswith('.gz'):
        with gzip.open(filepath, 'rt') as f:
            return json.load(f)
    else:
        with open(filepath, 'r') as f:
            return json.load(f)
