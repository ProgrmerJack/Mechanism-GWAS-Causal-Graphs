"""
Validate External Data Downloads
==================================

This script checks that all required external datasets for NBT validation
have been properly downloaded and are in the expected format.

Usage:
    python scripts/validate_external_data.py

Author: Mechanism-GWAS-Causal-Graphs team
Date: December 2025
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple
import json

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


class DataValidator:
    """Validate external dataset downloads."""
    
    def __init__(self, data_dir: Path = Path("data/external")):
        self.data_dir = data_dir
        self.validation_results = {}
        
    def check_file_exists(self, filepath: Path) -> Tuple[bool, str]:
        """Check if file exists and return status message."""
        if filepath.exists():
            size_mb = filepath.stat().st_size / (1024 * 1024)
            return True, f"✓ Found ({size_mb:.1f} MB)"
        else:
            return False, "✗ Missing"
    
    def validate_pops_data(self) -> Dict[str, str]:
        """Validate PoPS repository and data files."""
        pops_dir = self.data_dir / "pops"
        
        results = {
            "PoPS repository": self.check_file_exists(pops_dir),
            "Gene features": self.check_file_exists(pops_dir / "data" / "gene_features.txt"),
            "MAGMA scores": self.check_file_exists(pops_dir / "data" / "magma_0.05.genes.raw"),
        }
        
        return results
    
    def validate_crispr_data(self) -> Dict[str, str]:
        """Validate CRISPR validation datasets."""
        crispr_dir = self.data_dir / "crispr_validation"
        
        # Check for Fulco 2019 data (flexible naming)
        fulco_patterns = [
            "fulco_2019_table_s6a.xlsx",
            "41588_2019_538_MOESM3_ESM.xlsx",
            "fulco_validation.xlsx",
            "abc_validation.xlsx"
        ]
        
        fulco_found = False
        fulco_msg = "✗ Missing"
        for pattern in fulco_patterns:
            filepath = crispr_dir / pattern
            if filepath.exists():
                fulco_found = True
                size_mb = filepath.stat().st_size / (1024 * 1024)
                fulco_msg = f"✓ Found as {pattern} ({size_mb:.1f} MB)"
                break
        
        # Check for Gasperini 2019 data
        gasperini_patterns = [
            "gasperini_2019_screen.rds",
            "GSE120861_at_scale_screen.cds.rds",
            "gasperini_validation.rds"
        ]
        
        gasperini_found = False
        gasperini_msg = "✗ Missing"
        for pattern in gasperini_patterns:
            filepath = crispr_dir / pattern
            if filepath.exists():
                gasperini_found = True
                size_mb = filepath.stat().st_size / (1024 * 1024)
                gasperini_msg = f"✓ Found as {pattern} ({size_mb:.1f} MB)"
                break
        
        results = {
            "Fulco 2019 CRISPRi": (fulco_found, fulco_msg),
            "Gasperini 2019 screen": (gasperini_found, gasperini_msg),
        }
        
        return results
    
    def validate_drug_target_data(self) -> Dict[str, str]:
        """Validate drug-target datasets."""
        drug_dir = self.data_dir / "drug_targets"
        
        # Check ChEMBL (flexible formats)
        chembl_patterns = [
            "chembl_34.db",
            "chembl_33.db",
            "chembl.db",
            "approved_targets.tsv",
            "chembl_approved_targets.csv"
        ]
        
        chembl_found = False
        chembl_msg = "✗ Missing"
        for pattern in chembl_patterns:
            filepath = drug_dir / pattern
            if filepath.exists():
                chembl_found = True
                size_mb = filepath.stat().st_size / (1024 * 1024)
                chembl_msg = f"✓ Found as {pattern} ({size_mb:.1f} MB)"
                break
        
        # Check Open Targets
        ot_patterns = [
            "known_drugs.jsonl",
            "opentargets_drugs.json",
            "drug_targets.tsv"
        ]
        
        ot_found = False
        ot_msg = "✗ Missing (optional)"
        for pattern in ot_patterns:
            filepath = drug_dir / pattern
            if filepath.exists():
                ot_found = True
                size_mb = filepath.stat().st_size / (1024 * 1024)
                ot_msg = f"✓ Found as {pattern} ({size_mb:.1f} MB)"
                break
        
        results = {
            "ChEMBL drug targets": (chembl_found, chembl_msg),
            "Open Targets drugs": (ot_found, ot_msg),
        }
        
        return results
    
    def validate_open_targets_l2g(self) -> Dict[str, str]:
        """Validate Open Targets L2G scores."""
        ot_dir = self.data_dir / "open_targets"
        
        l2g_patterns = [
            "l2g_scores.parquet",
            "l2g_scores.tsv",
            "opentargets_l2g.parquet"
        ]
        
        l2g_found = False
        l2g_msg = "✗ Missing"
        for pattern in l2g_patterns:
            filepath = ot_dir / pattern
            if filepath.exists():
                l2g_found = True
                size_mb = filepath.stat().st_size / (1024 * 1024)
                l2g_msg = f"✓ Found as {pattern} ({size_mb:.1f} MB)"
                break
        
        results = {
            "Open Targets L2G scores": (l2g_found, l2g_msg),
        }
        
        return results
    
    def validate_flames_data(self) -> Dict[str, str]:
        """Validate FLAMES training data (optional)."""
        flames_dir = self.data_dir / "flames"
        
        flames_patterns = [
            "training_loci.tsv",
            "flames_training.csv",
            "exwas_loci.tsv"
        ]
        
        flames_found = False
        flames_msg = "✗ Missing (optional - will use approximation)"
        for pattern in flames_patterns:
            filepath = flames_dir / pattern
            if filepath.exists():
                flames_found = True
                size_mb = filepath.stat().st_size / (1024 * 1024)
                flames_msg = f"✓ Found as {pattern} ({size_mb:.1f} MB)"
                break
        
        results = {
            "FLAMES training data": (flames_found, flames_msg),
        }
        
        return results
    
    def run_all_validations(self) -> bool:
        """Run all validations and print report."""
        print("=" * 80)
        print("External Data Validation Report")
        print("=" * 80)
        print()
        
        all_valid = True
        
        # Critical datasets
        print("CRITICAL DATASETS (Required for NBT submission)")
        print("-" * 80)
        
        print("\n1. PoPS Baseline Method:")
        pops_results = self.validate_pops_data()
        for name, (valid, msg) in pops_results.items():
            print(f"   {name:30s} {msg}")
            if not valid and "optional" not in msg.lower():
                all_valid = False
        
        print("\n2. CRISPR Functional Validation:")
        crispr_results = self.validate_crispr_data()
        for name, (valid, msg) in crispr_results.items():
            print(f"   {name:30s} {msg}")
            if not valid and "optional" not in msg.lower():
                all_valid = False
        
        print("\n3. Drug-Target Validation:")
        drug_results = self.validate_drug_target_data()
        for name, (valid, msg) in drug_results.items():
            print(f"   {name:30s} {msg}")
            if not valid and "optional" not in msg.lower():
                all_valid = False
        
        print("\n4. Open Targets L2G Baseline:")
        l2g_results = self.validate_open_targets_l2g()
        for name, (valid, msg) in l2g_results.items():
            print(f"   {name:30s} {msg}")
            if not valid and "optional" not in msg.lower():
                all_valid = False
        
        # Optional datasets
        print("\n" + "=" * 80)
        print("OPTIONAL DATASETS (Nice to have)")
        print("-" * 80)
        
        print("\n5. FLAMES Method Data:")
        flames_results = self.validate_flames_data()
        for name, (valid, msg) in flames_results.items():
            print(f"   {name:30s} {msg}")
        
        # Summary
        print("\n" + "=" * 80)
        print("VALIDATION SUMMARY")
        print("=" * 80)
        
        if all_valid:
            print("✓ All critical datasets found!")
            print("  You are ready to proceed with baseline implementation.")
            print("\n  Next steps:")
            print("  1. Run: python src/baselines/baseline_runner.py")
            print("  2. See: NBT_UPGRADE_STRATEGY.md for full roadmap")
        else:
            print("✗ Some critical datasets are missing.")
            print("\n  Required actions:")
            print("  1. Review EXTERNAL_DATA_GUIDE.md for download instructions")
            print("  2. Download missing datasets")
            print("  3. Re-run this validation script")
            print("\n  Missing data will prevent NBT-grade validation.")
        
        print("=" * 80)
        
        return all_valid
    
    def generate_data_manifest(self, output_file: Path = Path("data/MANIFEST.json")):
        """Generate JSON manifest of available data."""
        manifest = {
            "validation_date": str(Path(__file__).stat().st_mtime),
            "datasets": {
                "pops": self.validate_pops_data(),
                "crispr": self.validate_crispr_data(),
                "drug_targets": self.validate_drug_target_data(),
                "open_targets_l2g": self.validate_open_targets_l2g(),
                "flames": self.validate_flames_data(),
            }
        }
        
        # Convert to JSON-serializable format
        serializable_manifest = {
            "validation_date": manifest["validation_date"],
            "datasets": {
                category: {
                    name: {"valid": valid, "message": msg}
                    for name, (valid, msg) in datasets.items()
                }
                for category, datasets in manifest["datasets"].items()
            }
        }
        
        with open(output_file, 'w') as f:
            json.dump(serializable_manifest, f, indent=2)
        
        print(f"\n✓ Data manifest written to: {output_file}")


def main():
    """Main validation entry point."""
    # Ensure data directory exists
    data_dir = Path("data/external")
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories if they don't exist
    for subdir in ["pops", "crispr_validation", "drug_targets", "open_targets", "flames"]:
        (data_dir / subdir).mkdir(exist_ok=True)
    
    # Run validation
    validator = DataValidator(data_dir)
    all_valid = validator.run_all_validations()
    
    # Generate manifest
    validator.generate_data_manifest()
    
    # Exit code
    sys.exit(0 if all_valid else 1)


if __name__ == "__main__":
    main()
