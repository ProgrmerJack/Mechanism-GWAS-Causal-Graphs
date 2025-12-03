"""
Logging utilities.
"""

import sys
import logging
from pathlib import Path
from typing import Optional
from datetime import datetime


def setup_logger(
    name: str = "mechanism_gwas",
    log_file: Optional[str] = None,
    level: str = "INFO",
    format_str: Optional[str] = None,
) -> logging.Logger:
    """
    Set up a logger with console and optional file handlers.
    
    Parameters
    ----------
    name : str
        Logger name.
    log_file : str, optional
        Path to log file.
    level : str
        Logging level (DEBUG, INFO, WARNING, ERROR).
    format_str : str, optional
        Log message format.
        
    Returns
    -------
    logging.Logger
        Configured logger.
    """
    if format_str is None:
        format_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    
    # Create logger
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper()))
    
    # Clear existing handlers
    logger.handlers = []
    
    # Create formatter
    formatter = logging.Formatter(format_str)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler (optional)
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        
        file_handler = logging.FileHandler(log_path)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger


def get_logger(name: str = "mechanism_gwas") -> logging.Logger:
    """
    Get an existing logger by name.
    
    Parameters
    ----------
    name : str
        Logger name.
        
    Returns
    -------
    logging.Logger
        Logger instance.
    """
    return logging.getLogger(name)


class ProgressLogger:
    """
    Simple progress logger for long-running tasks.
    """
    
    def __init__(
        self,
        total: int,
        desc: str = "Processing",
        logger: Optional[logging.Logger] = None,
        log_every: int = 100,
    ):
        """
        Initialize progress logger.
        
        Parameters
        ----------
        total : int
            Total number of items.
        desc : str
            Description of the task.
        logger : logging.Logger, optional
            Logger to use.
        log_every : int
            Log progress every N items.
        """
        self.total = total
        self.desc = desc
        self.logger = logger or get_logger()
        self.log_every = log_every
        self.current = 0
        self.start_time = datetime.now()
    
    def update(self, n: int = 1):
        """Update progress by n items."""
        self.current += n
        
        if self.current % self.log_every == 0 or self.current == self.total:
            elapsed = (datetime.now() - self.start_time).total_seconds()
            rate = self.current / elapsed if elapsed > 0 else 0
            pct = 100 * self.current / self.total
            
            self.logger.info(
                f"{self.desc}: {self.current}/{self.total} ({pct:.1f}%) - "
                f"{rate:.1f} items/sec"
            )
    
    def close(self):
        """Close progress logger and print summary."""
        elapsed = (datetime.now() - self.start_time).total_seconds()
        self.logger.info(
            f"{self.desc} complete: {self.current} items in {elapsed:.1f}s"
        )
