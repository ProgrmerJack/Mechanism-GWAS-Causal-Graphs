#!/usr/bin/env python3
"""
Data Download Script for Mechanism-GWAS-Causal-Graphs

Downloads all required datasets from public repositories:
- GWAS Catalog summary statistics
- ENCODE cCRE annotations
- GTEx v8 QTL summaries
- eQTL Catalogue data
- 1000 Genomes LD reference panels

Usage:
    python scripts/download_data.py --config config/config.yaml
    python scripts/download_data.py --dataset gwas_catalog
    python scripts/download_data.py --dataset encode_ccre
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
import yaml
from tqdm import tqdm


# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.utils.config import load_config
from src.utils.logging import setup_logger


logger = setup_logger("download_data")


class DataDownloader:
    """
    Handles downloading of all required datasets.
    """
    
    def __init__(self, config_path: str):
        """
        Initialize downloader with configuration.
        
        Parameters
        ----------
        config_path : str
            Path to configuration file.
        """
        self.config = load_config(config_path)
        self.project_root = Path(__file__).parent.parent
        self.data_dir = self.project_root / self.config["directories"]["data_raw"]
        
        # Create directories
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
    def download_file(
        self,
        url: str,
        output_path: Path,
        description: str = "",
        chunk_size: int = 8192,
    ) -> bool:
        """
        Download a file with progress bar.
        
        Parameters
        ----------
        url : str
            URL to download.
        output_path : Path
            Output file path.
        description : str
            Description for progress bar.
        chunk_size : int
            Download chunk size.
            
        Returns
        -------
        bool
            True if successful.
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Skip if already exists
        if output_path.exists():
            logger.info(f"File already exists: {output_path}")
            return True
        
        try:
            response = requests.get(url, stream=True, timeout=30)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            
            with open(output_path, 'wb') as f:
                with tqdm(
                    total=total_size,
                    unit='iB',
                    unit_scale=True,
                    desc=description or output_path.name,
                ) as pbar:
                    for chunk in response.iter_content(chunk_size=chunk_size):
                        size = f.write(chunk)
                        pbar.update(size)
            
            logger.info(f"Downloaded: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to download {url}: {e}")
            if output_path.exists():
                output_path.unlink()
            return False
    
    def download_gwas_catalog(self) -> None:
        """
        Download GWAS Catalog summary statistics for configured traits.
        """
        logger.info("=" * 60)
        logger.info("Downloading GWAS Catalog Summary Statistics")
        logger.info("=" * 60)
        
        output_dir = self.data_dir / "gwas_catalog"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        base_url = self.config["resources"]["gwas_catalog"]["harmonised_url"]
        
        for trait_config in self.config["gwas"]["traits"]:
            study_id = trait_config["study_id"]
            trait_name = trait_config["name"]
            
            # Construct URL (GWAS Catalog format)
            # Format: GCST90002232/harmonised/GCST90002232.h.tsv.gz
            url = f"{base_url}/{study_id}/harmonised/{study_id}.h.tsv.gz"
            output_file = output_dir / f"{study_id}_{trait_name}.h.tsv.gz"
            
            logger.info(f"Downloading {trait_name} ({study_id})")
            
            if not self.download_file(url, output_file, trait_name):
                # Try alternative URL format
                alt_url = f"{base_url}/{study_id}/{study_id}.h.tsv.gz"
                self.download_file(alt_url, output_file, trait_name)
        
        # Download metadata
        metadata_url = "https://www.ebi.ac.uk/gwas/api/search/downloads/full"
        metadata_file = output_dir / "gwas_catalog_metadata.tsv"
        self.download_file(metadata_url, metadata_file, "GWAS Catalog Metadata")
    
    def download_encode_ccre(self) -> None:
        """
        Download ENCODE cCRE annotations from SCREEN.
        """
        logger.info("=" * 60)
        logger.info("Downloading ENCODE cCRE Annotations")
        logger.info("=" * 60)
        
        output_dir = self.data_dir / "encode_ccre"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # SCREEN cCRE files
        ccre_files = {
            # Human cCREs (GRCh38)
            "GRCh38-cCREs.bed": "https://api.wenglab.org/screen_v13/fdownloads/GRCh38-cCREs.bed",
            
            # Biosample-specific activity
            "GRCh38-cCREs.CTCF-only.bed": "https://api.wenglab.org/screen_v13/fdownloads/cCREs/GRCh38-CTCF-only.bed",
            "GRCh38-cCREs.dELS.bed": "https://api.wenglab.org/screen_v13/fdownloads/cCREs/GRCh38-dELS.bed",
            "GRCh38-cCREs.pELS.bed": "https://api.wenglab.org/screen_v13/fdownloads/cCREs/GRCh38-pELS.bed",
            "GRCh38-cCREs.PLS.bed": "https://api.wenglab.org/screen_v13/fdownloads/cCREs/GRCh38-PLS.bed",
        }
        
        for filename, url in ccre_files.items():
            output_file = output_dir / filename
            self.download_file(url, output_file, filename)
        
        # Download tissue/biosample-specific activity matrices if available
        # These are larger files, download selectively
        logger.info("Note: For biosample-specific cCRE activity, use ENCODE Portal or AWS")
        
        # Create a helper script for AWS download
        aws_script = output_dir / "download_encode_aws.sh"
        with open(aws_script, 'w') as f:
            f.write("""#!/bin/bash
# Download ENCODE data from AWS (faster for large files)
# Requires: aws cli configured with anonymous access

# cCRE files
aws s3 cp --no-sign-request \\
    s3://encode-public/2022/11/14/de0f8e8e-6fb7-4720-b3b0-8e95e1e6ec8c/GRCh38-cCREs.bed \\
    GRCh38-cCREs.bed

# Add more files as needed
""")
        
        logger.info(f"AWS download script created: {aws_script}")
    
    def download_gtex(self) -> None:
        """
        Download GTEx v8 QTL summaries.
        """
        logger.info("=" * 60)
        logger.info("Downloading GTEx v8 QTL Data")
        logger.info("=" * 60)
        
        output_dir = self.data_dir / "gtex_v8"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        base_url = self.config["resources"]["gtex"]["ftp_url"]
        
        # Download eQTL files for configured tissues
        tissues = self.config["qtl"]["gtex"]["tissues"]
        
        for tissue in tissues:
            # eQTL significant pairs
            eqtl_url = f"{base_url}/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL/{tissue}.v8.signif_variant_gene_pairs.txt.gz"
            eqtl_file = output_dir / f"{tissue}.v8.signif_variant_gene_pairs.txt.gz"
            
            self.download_file(eqtl_url, eqtl_file, f"GTEx eQTL - {tissue}")
            
            # sQTL significant pairs
            sqtl_url = f"{base_url}/single-tissue-cis-qtl/GTEx_Analysis_v8_sQTL/{tissue}.v8.sqtl_signifpairs.txt.gz"
            sqtl_file = output_dir / f"{tissue}.v8.sqtl_signifpairs.txt.gz"
            
            self.download_file(sqtl_url, sqtl_file, f"GTEx sQTL - {tissue}")
        
        # Download gene annotation
        gene_url = f"{base_url}/reference/gencode.v26.GRCh38.genes.gtf"
        gene_file = output_dir / "gencode.v26.GRCh38.genes.gtf"
        self.download_file(gene_url, gene_file, "GTEx Gene Annotation")
        
        # Download lookup tables
        lookup_url = f"{base_url}/reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"
        lookup_file = output_dir / "GTEx_lookup_table.txt.gz"
        self.download_file(lookup_url, lookup_file, "GTEx Lookup Table")
    
    def download_eqtl_catalogue(self) -> None:
        """
        Download eQTL Catalogue data.
        """
        logger.info("=" * 60)
        logger.info("Downloading eQTL Catalogue Data")
        logger.info("=" * 60)
        
        output_dir = self.data_dir / "eqtl_catalogue"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Download study metadata
        api_url = self.config["resources"]["eqtl_catalogue"]["api_url"]
        
        # Get list of available datasets
        try:
            response = requests.get(f"{api_url}/datasets", timeout=30)
            response.raise_for_status()
            datasets = response.json()
            
            # Save dataset list
            with open(output_dir / "datasets.json", 'w') as f:
                import json
                json.dump(datasets, f, indent=2)
            
            logger.info(f"Found {len(datasets)} datasets in eQTL Catalogue")
            
        except Exception as e:
            logger.error(f"Failed to fetch eQTL Catalogue metadata: {e}")
        
        # Download fine-mapped results (these are the most useful)
        ftp_url = self.config["resources"]["eqtl_catalogue"]["ftp_url"]
        
        # Priority studies (large, high-quality)
        priority_studies = [
            "GTEx_V8",  # Will be separate from our GTEx download
            "BLUEPRINT",
            "TwinsUK",
            "Alasoo_2018",
        ]
        
        # Create download manifest
        manifest_content = f"""# eQTL Catalogue Download Manifest
# FTP Base: {ftp_url}

# Fine-mapped results (SuSiE credible sets)
# Format: ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/[dataset_id]/

# Download example:
# wget -r -np -nH --cut-dirs=5 {ftp_url}/susie/GTEx_V8/

# Note: Full catalogue is ~500GB. Download selectively based on tissue relevance.
"""
        
        with open(output_dir / "DOWNLOAD_MANIFEST.md", 'w') as f:
            f.write(manifest_content)
        
        logger.info("eQTL Catalogue manifest created. Download specific studies as needed.")
    
    def download_ld_reference(self) -> None:
        """
        Download 1000 Genomes LD reference panels.
        """
        logger.info("=" * 60)
        logger.info("Downloading 1000 Genomes LD Reference Panels")
        logger.info("=" * 60)
        
        output_dir = self.data_dir / "ld_reference"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Use the PLINK-formatted 1000G files
        # These are available from various sources
        
        ld_url = self.config["resources"]["ld_reference"]["url"]
        
        # Create download script for large files
        download_script = output_dir / "download_ld_panels.sh"
        
        script_content = """#!/bin/bash
# Download 1000 Genomes Phase 3 LD Reference Panels
# These files are large (~15GB per population)

# Base URL
BASE_URL="https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b38.vcf"

# Populations
POPS=("EUR" "EAS" "AFR" "AMR" "SAS")

# Download per-chromosome VCF files
for pop in "${POPS[@]}"; do
    mkdir -p ${pop}
    for chr in $(seq 1 22); do
        echo "Downloading ${pop} chr${chr}..."
        wget -c -P ${pop}/ "${BASE_URL}/chr${chr}.1kg.phase3.v5a.b38.vcf.gz"
        wget -c -P ${pop}/ "${BASE_URL}/chr${chr}.1kg.phase3.v5a.b38.vcf.gz.tbi"
    done
done

echo "Download complete!"
echo "Next: Convert to PLINK format with:"
echo "  plink2 --vcf chr1.vcf.gz --make-bed --out chr1"
"""
        
        with open(download_script, 'w') as f:
            f.write(script_content)
        
        logger.info(f"LD reference download script created: {download_script}")
        logger.info("Run the script manually due to large file sizes (~100GB total)")
        
        # Download sample information
        sample_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index"
        sample_file = output_dir / "1000G_sample_info.tsv"
        
        # Download population panel
        panel_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_all_populations.vcf.gz"
        
        logger.info("Note: For full LD panels, run the download script in the ld_reference directory")
    
    def download_benchmark_data(self) -> None:
        """
        Download benchmark/validation datasets.
        """
        logger.info("=" * 60)
        logger.info("Downloading Benchmark Datasets")
        logger.info("=" * 60)
        
        output_dir = self.data_dir.parent / "external"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # OMIM - Mendelian disease genes (requires registration)
        omim_readme = output_dir / "OMIM_README.md"
        with open(omim_readme, 'w') as f:
            f.write("""# OMIM Data

To download OMIM gene-disease associations:
1. Register at https://omim.org/downloads
2. Download genemap2.txt
3. Place in this directory

Alternative: Use the OMIM API with your API key.
""")
        
        # DrugBank targets (requires registration)
        drugbank_readme = output_dir / "DrugBank_README.md"
        with open(drugbank_readme, 'w') as f:
            f.write("""# DrugBank Data

To download DrugBank drug-target associations:
1. Register at https://go.drugbank.com/
2. Download the approved drug targets file
3. Place in this directory

Used for: Benchmarking gene prioritization against known drug targets
""")
        
        # Open Targets L2G scores (for comparison)
        ot_url = "https://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/lut/variant-index/"
        ot_readme = output_dir / "OpenTargets_README.md"
        with open(ot_readme, 'w') as f:
            f.write(f"""# Open Targets Genetics Data

Download Locus-to-Gene (L2G) scores for comparison:
{ot_url}

These are our primary baseline for benchmarking.
""")
        
        logger.info("Benchmark data readme files created")
        logger.info("Manual registration required for some resources (OMIM, DrugBank)")
    
    def download_all(self) -> None:
        """
        Download all required datasets.
        """
        logger.info("Starting full data download...")
        
        self.download_gwas_catalog()
        self.download_encode_ccre()
        self.download_gtex()
        self.download_eqtl_catalogue()
        self.download_ld_reference()
        self.download_benchmark_data()
        
        logger.info("=" * 60)
        logger.info("Data download complete!")
        logger.info("=" * 60)
        logger.info(f"Data directory: {self.data_dir}")
        logger.info("Note: Some large files require manual download (see scripts)")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Download data for Mechanism-GWAS-Causal-Graphs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    parser.add_argument(
        "--config",
        type=str,
        default="config/config.yaml",
        help="Path to configuration file",
    )
    
    parser.add_argument(
        "--dataset",
        type=str,
        choices=["gwas_catalog", "encode_ccre", "gtex", "eqtl_catalogue", "ld_reference", "benchmark", "all"],
        default="all",
        help="Specific dataset to download",
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print URLs without downloading",
    )
    
    args = parser.parse_args()
    
    # Find config file
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = Path(__file__).parent.parent / config_path
    
    if not config_path.exists():
        logger.error(f"Configuration file not found: {config_path}")
        sys.exit(1)
    
    # Initialize downloader
    downloader = DataDownloader(str(config_path))
    
    # Download requested dataset(s)
    if args.dataset == "all":
        downloader.download_all()
    elif args.dataset == "gwas_catalog":
        downloader.download_gwas_catalog()
    elif args.dataset == "encode_ccre":
        downloader.download_encode_ccre()
    elif args.dataset == "gtex":
        downloader.download_gtex()
    elif args.dataset == "eqtl_catalogue":
        downloader.download_eqtl_catalogue()
    elif args.dataset == "ld_reference":
        downloader.download_ld_reference()
    elif args.dataset == "benchmark":
        downloader.download_benchmark_data()


if __name__ == "__main__":
    main()
