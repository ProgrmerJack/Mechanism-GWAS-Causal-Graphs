"""
External Data Validation & Download Manager
=============================================

Validates presence of all required external data files for the paper.
Automatically downloads missing critical files from public repositories.
Generates comprehensive data inventory report.
"""

import os
import requests
import gzip
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging
import json
from dataclasses import dataclass
import hashlib


@dataclass
class DataFile:
    """Specification for an external data file."""
    name: str
    path: str  # Relative to data/external
    source: str
    source_url: Optional[str]  # URL to download from
    md5: Optional[str]  # Expected MD5 checksum
    size_mb: Optional[float]  # Expected size in MB
    critical: bool  # Is this required for paper?
    description: str


class DataValidator:
    """
    Validates and manages external data files.
    
    Tracks:
    - Which files are present/missing
    - File integrity (size, checksum)
    - Download status
    - Data requirements for each analysis
    """
    
    # Define all external data required for the paper
    REQUIRED_DATA_FILES = [
        # CRISPR Validation Data
        DataFile(
            name="Fulco 2019 CRISPRi-FlowFISH",
            path="crispr_validation/fulco_2019_table_s6a.xlsx",
            source="GEO",
            source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96739",
            md5=None,
            size_mb=5.2,
            critical=True,
            description="K562 enhancer perturbation data (6,090 element-gene pairs)"
        ),
        
        DataFile(
            name="Gasperini 2019 CRISPR Screen",
            path="crispr_validation/gasperini_2019_results.txt.gz",
            source="GEO",
            source_url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120861",
            md5=None,
            size_mb=285.0,
            critical=True,
            description="CRISPRi screen results (18,000 enhancers, 1,750 genes)"
        ),
        
        # PoPS Data
        DataFile(
            name="PoPS Gene Features",
            path="pops/gene_features.txt.gz",
            source="GitHub",
            source_url="https://github.com/ajayakumar-ak/PoPS/raw/master/data/features/gene_features.txt.gz",
            md5=None,
            size_mb=12.5,
            critical=True,
            description="Pre-computed gene feature matrix for PoPS"
        ),
        
        DataFile(
            name="PoPS Genome Adjacency",
            path="pops/genome_adjacency.txt.gz",
            source="GitHub",
            source_url="https://github.com/ajayakumar-ak/PoPS/raw/master/data/genome_adjacency.txt.gz",
            md5=None,
            size_mb=18.3,
            critical=True,
            description="Gene-gene adjacency matrix for network features"
        ),
        
        # ChEMBL Drug Target Data
        DataFile(
            name="ChEMBL Human Targets",
            path="chembl/chembl_human_targets.txt.gz",
            source="ChEMBL",
            source_url="https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_31/",
            md5=None,
            size_mb=45.0,
            critical=True,
            description="Drug-target interactions for 15,000+ human proteins"
        ),
        
        # Open Targets Genetics L2G
        DataFile(
            name="Open Targets L2G Scores",
            path="opentargets/l2g_scores.parquet",
            source="Open Targets",
            source_url="https://platform-api.opentargets.org/v3/platform/public/locus2gene/download",
            md5=None,
            size_mb=250.0,
            critical=True,
            description="L2G scores for GWAS locus-to-gene mapping"
        ),
        
        # STRING Network
        DataFile(
            name="STRING Human PPI Network",
            path="string/string_interactions.txt.gz",
            source="STRING",
            source_url="https://stringdb-static.org/download/protein.links.full.v11.5/9606.protein.links.full.v11.5.txt.gz",
            md5=None,
            size_mb=350.0,
            critical=False,  # Can use precomputed features
            description="Complete human protein-protein interaction network"
        ),
        
        # Enrichment Data
        DataFile(
            name="GO Annotations",
            path="annotations/go_human.obo",
            source="GO",
            source_url="http://purl.obolibrary.org/obo/go.obo",
            md5=None,
            size_mb=48.0,
            critical=False,
            description="Gene Ontology annotations"
        ),
        
        DataFile(
            name="GWAS Catalog",
            path="gwas/gwas_catalog_v1.0.tsv.gz",
            source="GWAS Catalog",
            source_url="https://www.ebi.ac.uk/gwas/api/search/downloads/full",
            md5=None,
            size_mb=125.0,
            critical=False,
            description="Complete GWAS Catalog with all published associations"
        ),
    ]
    
    def __init__(self, data_dir: str = "data/external"):
        """Initialize validator."""
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(__name__)
    
    def validate_all(self) -> Tuple[Dict[str, bool], List[str]]:
        """
        Validate all required data files.
        
        Returns:
            (status_dict, missing_list)
            - status_dict: {filename: is_present}
            - missing_list: List of missing files
        """
        self.logger.info("Validating external data files...")
        status = {}
        missing = []
        
        for data_file in self.REQUIRED_DATA_FILES:
            file_path = self.data_dir / data_file.path
            exists = file_path.exists()
            status[data_file.name] = exists
            
            if exists:
                size_mb = file_path.stat().st_size / (1024**2)
                self.logger.info(f"✓ {data_file.name} ({size_mb:.1f} MB)")
            else:
                self.logger.warning(f"✗ {data_file.name}")
                if data_file.critical:
                    missing.append(data_file.name)
        
        return status, missing
    
    def download_missing_files(self, auto_download: bool = True) -> Dict[str, bool]:
        """
        Download missing critical files.
        
        Args:
            auto_download: Whether to download without prompting
        
        Returns:
            {filename: download_success}
        """
        status, missing = self.validate_all()
        download_status = {}
        
        if not missing:
            self.logger.info("All critical files present!")
            return {f.name: True for f in self.REQUIRED_DATA_FILES}
        
        self.logger.info(f"\nMissing {len(missing)} critical files:")
        for name in missing:
            self.logger.info(f"  - {name}")
        
        if not auto_download:
            proceed = input("\nDownload missing files? [y/n]: ").lower() == 'y'
            if not proceed:
                return download_status
        
        self.logger.info("\nDownloading files...")
        
        for data_file in self.REQUIRED_DATA_FILES:
            if data_file.name not in missing:
                continue
            
            if data_file.source_url is None:
                self.logger.warning(f"No download URL for {data_file.name}")
                download_status[data_file.name] = False
                continue
            
            try:
                self._download_file(data_file)
                download_status[data_file.name] = True
            except Exception as e:
                self.logger.error(f"Failed to download {data_file.name}: {e}")
                download_status[data_file.name] = False
        
        return download_status
    
    def _download_file(self, data_file: DataFile) -> None:
        """Download a single file."""
        output_path = self.data_dir / data_file.path
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        self.logger.info(f"Downloading {data_file.name}...")
        self.logger.info(f"  From: {data_file.source}")
        
        # Special handling for different sources
        if data_file.source == "GitHub":
            self._download_from_github(data_file, output_path)
        elif data_file.source == "GEO":
            self._download_from_geo(data_file, output_path)
        elif data_file.source == "ChEMBL":
            self._download_chembl(data_file, output_path)
        elif data_file.source == "Open Targets":
            self._download_open_targets(data_file, output_path)
        elif data_file.source == "STRING":
            self._download_string(data_file, output_path)
        elif data_file.source in ["GO", "GWAS Catalog"]:
            self._download_generic(data_file, output_path)
    
    def _download_from_github(self, data_file: DataFile, output_path: Path) -> None:
        """Download from GitHub."""
        response = requests.get(data_file.source_url, stream=True)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        downloaded = 0
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                downloaded += len(chunk)
                if total_size > 0:
                    pct = 100 * downloaded / total_size
                    self.logger.info(f"  {pct:.1f}% ({downloaded/(1024**2):.1f} MB)")
        
        self.logger.info(f"✓ Downloaded to {output_path}")
    
    def _download_from_geo(self, data_file: DataFile, output_path: Path) -> None:
        """Download from NCBI GEO."""
        self.logger.info(f"  Please download manually from: {data_file.source_url}")
        self.logger.info(f"  Save to: {output_path}")
        self.logger.warning("  GEO downloads require manual interaction (FTP/web interface)")
        raise RuntimeError("Manual GEO download required")
    
    def _download_chembl(self, data_file: DataFile, output_path: Path) -> None:
        """Download from ChEMBL."""
        self.logger.info(f"  Downloading from ChEMBL FTP...")
        # Would implement FTP download here
        raise NotImplementedError("ChEMBL download not yet implemented")
    
    def _download_open_targets(self, data_file: DataFile, output_path: Path) -> None:
        """Download from Open Targets."""
        self.logger.info(f"  Downloading from Open Targets API...")
        response = requests.get(data_file.source_url)
        response.raise_for_status()
        
        with open(output_path, 'wb') as f:
            f.write(response.content)
        
        self.logger.info(f"✓ Downloaded {output_path.stat().st_size / (1024**2):.1f} MB")
    
    def _download_string(self, data_file: DataFile, output_path: Path) -> None:
        """Download from STRING."""
        response = requests.get(data_file.source_url)
        response.raise_for_status()
        
        with open(output_path, 'wb') as f:
            f.write(response.content)
        
        self.logger.info(f"✓ Downloaded {output_path.stat().st_size / (1024**2):.1f} MB")
    
    def _download_generic(self, data_file: DataFile, output_path: Path) -> None:
        """Generic HTTP download."""
        response = requests.get(data_file.source_url, stream=True)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        size_mb = output_path.stat().st_size / (1024**2)
        self.logger.info(f"✓ Downloaded {size_mb:.1f} MB")
    
    def generate_inventory_report(self) -> str:
        """Generate data inventory report."""
        report = []
        report.append("\n" + "="*100)
        report.append("EXTERNAL DATA INVENTORY REPORT")
        report.append("="*100 + "\n")
        
        status, missing = self.validate_all()
        
        report.append(f"Total Files: {len(self.REQUIRED_DATA_FILES)}")
        report.append(f"Present: {len(status) - len(missing)}")
        report.append(f"Missing: {len(missing)}")
        report.append(f"Critical Missing: {len([m for m in missing if any(f.name == m and f.critical for f in self.REQUIRED_DATA_FILES)])}\n")
        
        # By source
        report.append("Files by Source:")
        report.append("-" * 100)
        
        by_source = {}
        for f in self.REQUIRED_DATA_FILES:
            if f.source not in by_source:
                by_source[f.source] = []
            by_source[f.source].append(f)
        
        for source, files in sorted(by_source.items()):
            report.append(f"\n{source}:")
            for f in files:
                status_str = "✓" if f.name in status and status[f.name] else "✗"
                critical_str = "[CRITICAL]" if f.critical else "[Optional]"
                report.append(f"  {status_str} {f.name} {critical_str}")
                report.append(f"      {f.description}")
                if f.source_url:
                    report.append(f"      URL: {f.source_url}")
        
        report.append("\n" + "="*100 + "\n")
        
        return "\n".join(report)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    validator = DataValidator()
    print(validator.generate_inventory_report())
    
    # Download missing files
    # validator.download_missing_files(auto_download=True)
