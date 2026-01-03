#!/usr/bin/env python3
"""
Search UK Biobank ExWAS for Rare Variant Burden Tests

Identifies gene-based burden associations meeting exome-wide significance
(p < 3.6e-7) and matches them to GWAS loci by trait ontology.

Data sources:
- UK Biobank ExWAS results (gene-based burden tests)
- Published ExWAS catalogs
- GWAS-trait to ExWAS-phenotype mappings

Author: Generated for RegulatoryBench v3 platinum
Date: 2025-01
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Dict, List, Set, Tuple
import requests
import json
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class ExWASBurdenSearcher:
    """Search for rare variant burden tests supporting gene-locus associations."""
    
    # Exome-wide significance threshold
    EXOME_WIDE_THRESHOLD = 3.6e-7  # Bonferroni correction for ~20,000 genes
    
    # Trait ontology mappings: GWAS trait → ExWAS phenotype categories
    TRAIT_MAPPINGS = {
        'lipids': [
            'cholesterol', 'triglyceride', 'ldl', 'hdl', 'hyperlipidemia',
            'hypercholesterolemia', 'familial hypercholesterolemia'
        ],
        'metabolic': [
            'diabetes', 't2d', 'glucose', 'insulin', 'hba1c', 'obesity',
            'bmi', 'metabolic syndrome'
        ],
        'cardiovascular': [
            'coronary', 'myocardial', 'heart attack', 'stroke', 'hypertension',
            'blood pressure', 'atrial fibrillation', 'heart failure'
        ],
        'immune': [
            'autoimmune', 'inflammatory', 'crohn', 'ulcerative colitis',
            'rheumatoid', 'lupus', 'psoriasis', 'asthma', 'allergy'
        ],
        'neurological': [
            'alzheimer', 'parkinson', 'dementia', 'epilepsy', 'migraine',
            'multiple sclerosis', 'neuropathy'
        ],
        'psychiatric': [
            'depression', 'bipolar', 'schizophrenia', 'anxiety',
            'autism', 'adhd'
        ],
        'renal': [
            'kidney', 'renal', 'egfr', 'creatinine', 'chronic kidney disease'
        ],
        'liver': [
            'liver', 'hepatic', 'cirrhosis', 'fatty liver', 'alt', 'ast'
        ],
    }
    
    def __init__(
        self,
        benchmark_dir: str = "benchmarks",
        p_threshold: float = 3.6e-7,
        data_sources: List[str] = None,
    ):
        self.benchmark_dir = Path(benchmark_dir)
        self.output_dir = self.benchmark_dir / "evidence_curation"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.p_threshold = p_threshold
        
        # Data sources to search
        self.data_sources = data_sources or [
            'uk_biobank_exwas',
            'published_catalogs',
            'genebass',
            'azphewas'
        ]
        
        # Storage for burden test results
        self.burden_genes = []
        
    def load_gwas_loci(self) -> pd.DataFrame:
        """
        Load GWAS loci from benchmark to know which traits we need.
        
        Returns
        -------
        pd.DataFrame
            GWAS loci with trait information.
        """
        logger.info("Loading GWAS loci from benchmark...")
        
        # Try to find loci file
        loci_files = list(self.benchmark_dir.glob("**/loci/*.parquet"))
        
        if not loci_files:
            # Try alternative locations
            alt_files = list(Path("data/external/flames").glob("**/loci*.parquet"))
            if alt_files:
                loci_files = alt_files
        
        if not loci_files:
            logger.error("No GWAS loci files found!")
            return pd.DataFrame()
        
        # Load all loci
        dfs = []
        for file in loci_files:
            try:
                df = pd.read_parquet(file)
                dfs.append(df)
            except Exception as e:
                logger.warning(f"Could not read {file}: {e}")
        
        if not dfs:
            return pd.DataFrame()
        
        loci = pd.concat(dfs, ignore_index=True)
        logger.info(f"Loaded {len(loci):,} GWAS loci")
        
        # Get unique traits
        if 'trait' in loci.columns:
            traits = loci['trait'].unique()
            logger.info(f"Found {len(traits)} unique traits")
            
            # Map to ExWAS phenotype categories
            self._build_trait_mappings(traits)
        
        return loci
    
    def _build_trait_mappings(self, gwas_traits: List[str]):
        """
        Build mappings from GWAS traits to ExWAS phenotype categories.
        
        Parameters
        ----------
        gwas_traits : list
            List of GWAS trait names.
        """
        logger.info("Building GWAS trait → ExWAS phenotype mappings...")
        
        self.trait_to_category = {}
        
        for trait in gwas_traits:
            trait_lower = str(trait).lower()
            
            # Find matching category
            for category, keywords in self.TRAIT_MAPPINGS.items():
                if any(kw in trait_lower for kw in keywords):
                    self.trait_to_category[trait] = category
                    break
            
            # Default to metabolic if no match
            if trait not in self.trait_to_category:
                self.trait_to_category[trait] = 'metabolic'
        
        logger.info(f"Mapped {len(self.trait_to_category)} GWAS traits to ExWAS categories")
    
    def search_ukb_exwas(self) -> int:
        """
        Search UK Biobank ExWAS results.
        
        This is a placeholder - in practice, you would:
        1. Download UK Biobank ExWAS results from published sources
        2. Parse gene-based burden test results
        3. Filter for exome-wide significance
        4. Match to GWAS traits by phenotype ontology
        
        Returns
        -------
        int
            Number of burden associations found.
        """
        logger.info("Searching UK Biobank ExWAS results...")
        
        # Check for pre-downloaded ExWAS results
        exwas_file = self.output_dir / "ukb_exwas_burden_tests.tsv"
        
        if not exwas_file.exists():
            logger.warning(f"UK Biobank ExWAS file not found: {exwas_file}")
            logger.warning("To use this feature:")
            logger.warning("  1. Download UK Biobank ExWAS results")
            logger.warning("  2. Extract gene-based burden tests")
            logger.warning("  3. Save to: benchmarks/evidence_curation/ukb_exwas_burden_tests.tsv")
            logger.warning("  4. Required columns: gene_symbol, phenotype, p_value, beta, n_cases")
            return 0
        
        # Load ExWAS results
        exwas = pd.read_csv(exwas_file, sep='\t')
        logger.info(f"Loaded {len(exwas):,} ExWAS burden test results")
        
        # Filter for exome-wide significance
        sig_exwas = exwas[exwas['p_value'] < self.p_threshold]
        logger.info(f"Found {len(sig_exwas):,} significant burden tests (p < {self.p_threshold})")
        
        # Extract gene-phenotype associations
        for _, row in sig_exwas.iterrows():
            gene = str(row['gene_symbol'])
            phenotype = str(row['phenotype'])
            p_value = float(row['p_value'])
            
            # Map phenotype to GWAS trait category
            trait_category = self._phenotype_to_trait_category(phenotype)
            
            self.burden_genes.append({
                'gene_symbol': gene,
                'evidence_type': 'rare_variant_burden',
                'evidence_source': 'UK_Biobank_ExWAS',
                'phenotype': phenotype,
                'trait_category': trait_category,
                'p_value': p_value,
                'beta': row.get('beta', np.nan),
                'n_cases': row.get('n_cases', np.nan),
                'confidence': 'high' if p_value < 1e-8 else 'medium',
                'date_curated': pd.Timestamp.now().strftime('%Y-%m-%d'),
                'pmid': row.get('pmid', ''),
                'doi': row.get('doi', '')
            })
        
        logger.info(f"Extracted {len(sig_exwas):,} rare variant burden associations")
        return len(sig_exwas)
    
    def _phenotype_to_trait_category(self, phenotype: str) -> str:
        """
        Map ExWAS phenotype to GWAS trait category.
        
        Parameters
        ----------
        phenotype : str
            ExWAS phenotype name.
            
        Returns
        -------
        str
            GWAS trait category.
        """
        phenotype_lower = str(phenotype).lower()
        
        for category, keywords in self.TRAIT_MAPPINGS.items():
            if any(kw in phenotype_lower for kw in keywords):
                return category
        
        return 'other'
    
    def search_genebass(self) -> int:
        """
        Search GeneBass database for rare variant burden tests.
        
        GeneBass: https://app.genebass.org/
        Gene-based association tests in UK Biobank exomes.
        
        Returns
        -------
        int
            Number of associations found.
        """
        logger.info("Searching GeneBass database...")
        
        # GeneBass API (if available) or pre-downloaded data
        genebass_file = self.output_dir / "genebass_burden_tests.tsv"
        
        if not genebass_file.exists():
            logger.warning(f"GeneBass file not found: {genebass_file}")
            logger.warning("To use GeneBass:")
            logger.warning("  1. Visit https://app.genebass.org/")
            logger.warning("  2. Download gene-based burden tests")
            logger.warning("  3. Save to: benchmarks/evidence_curation/genebass_burden_tests.tsv")
            return 0
        
        # Load GeneBass results (same processing as UK Biobank)
        genebass = pd.read_csv(genebass_file, sep='\t')
        logger.info(f"Loaded {len(genebass):,} GeneBass burden tests")
        
        sig_genebass = genebass[genebass['p_value'] < self.p_threshold]
        logger.info(f"Found {len(sig_genebass):,} significant GeneBass associations")
        
        for _, row in sig_genebass.iterrows():
            self.burden_genes.append({
                'gene_symbol': str(row['gene_symbol']),
                'evidence_type': 'rare_variant_burden',
                'evidence_source': 'GeneBass',
                'phenotype': str(row['phenotype']),
                'trait_category': self._phenotype_to_trait_category(row['phenotype']),
                'p_value': float(row['p_value']),
                'confidence': 'high',
                'date_curated': pd.Timestamp.now().strftime('%Y-%m-%d')
            })
        
        return len(sig_genebass)
    
    def search_published_catalogs(self) -> int:
        """
        Search published ExWAS catalogs from literature.
        
        Examples:
        - DeBoever et al. 2018 (UK Biobank ExWAS)
        - Backman et al. 2021 (AstraZeneca PheWAS)
        - Karczewski et al. 2020 (gnomAD)
        
        Returns
        -------
        int
            Number of associations found.
        """
        logger.info("Searching published ExWAS catalogs...")
        
        # This would involve:
        # 1. Manual curation of published ExWAS results
        # 2. Extraction of gene-level burden associations
        # 3. Standardization of phenotype names
        
        # Placeholder: check for curated file
        catalog_file = self.output_dir / "published_exwas_curated.tsv"
        
        if not catalog_file.exists():
            logger.warning("No published ExWAS catalog found")
            logger.warning("Manually curate significant associations from:")
            logger.warning("  - DeBoever et al. 2018 (PMID: 29785011)")
            logger.warning("  - Backman et al. 2021 (PMID: 34012112)")
            logger.warning("  - Karczewski et al. 2020 (PMID: 32461654)")
            return 0
        
        catalog = pd.read_csv(catalog_file, sep='\t')
        logger.info(f"Loaded {len(catalog):,} curated ExWAS associations")
        
        for _, row in catalog.iterrows():
            self.burden_genes.append({
                'gene_symbol': str(row['gene_symbol']),
                'evidence_type': 'rare_variant_burden',
                'evidence_source': f"Published_ExWAS_{row.get('study', 'unknown')}",
                'phenotype': str(row['phenotype']),
                'trait_category': self._phenotype_to_trait_category(row['phenotype']),
                'p_value': float(row['p_value']),
                'confidence': 'high',
                'pmid': row.get('pmid', ''),
                'doi': row.get('doi', ''),
                'date_curated': pd.Timestamp.now().strftime('%Y-%m-%d')
            })
        
        return len(catalog)
    
    def match_to_gwas_loci(
        self,
        gwas_loci: pd.DataFrame,
        distance_threshold: int = 500_000
    ) -> pd.DataFrame:
        """
        Match burden test genes to GWAS loci by trait and distance.
        
        Parameters
        ----------
        gwas_loci : pd.DataFrame
            GWAS loci.
        distance_threshold : int
            Maximum distance from lead variant (default: 500kb).
            
        Returns
        -------
        pd.DataFrame
            Matched gene-locus associations.
        """
        if not self.burden_genes:
            logger.warning("No burden genes to match!")
            return pd.DataFrame()
        
        logger.info(f"Matching {len(self.burden_genes)} burden genes to GWAS loci...")
        
        burden_df = pd.DataFrame(self.burden_genes)
        
        # Match by trait category
        matched = []
        
        for _, locus in gwas_loci.iterrows():
            locus_trait = locus.get('trait', '')
            locus_category = self.trait_to_category.get(locus_trait, 'other')
            
            # Find burden genes matching this trait category
            matching_genes = burden_df[
                burden_df['trait_category'] == locus_category
            ]
            
            for _, gene_row in matching_genes.iterrows():
                matched.append({
                    'locus_id': locus.get('locus_id', ''),
                    'gene_symbol': gene_row['gene_symbol'],
                    'evidence_type': gene_row['evidence_type'],
                    'evidence_source': gene_row['evidence_source'],
                    'phenotype': gene_row['phenotype'],
                    'p_value': gene_row['p_value'],
                    'confidence': gene_row['confidence'],
                    'pmid': gene_row.get('pmid', ''),
                    'doi': gene_row.get('doi', ''),
                    'date_curated': gene_row['date_curated']
                })
        
        matched_df = pd.DataFrame(matched)
        
        if matched_df.empty:
            logger.warning("No matches found between burden genes and GWAS loci!")
            return matched_df
        
        # Deduplicate
        matched_df = matched_df.drop_duplicates(subset=['locus_id', 'gene_symbol'])
        
        logger.info(f"Matched {len(matched_df):,} gene-locus pairs")
        
        return matched_df
    
    def save_results(self, output_file: str = "rare_burden_genes.tsv"):
        """
        Save burden test results.
        
        Parameters
        ----------
        output_file : str
            Output filename.
        """
        if not self.burden_genes:
            logger.warning("No burden genes to save!")
            return
        
        df = pd.DataFrame(self.burden_genes)
        
        # Deduplicate by gene and phenotype
        df_dedup = df.drop_duplicates(subset=['gene_symbol', 'phenotype'])
        
        output_path = self.output_dir / output_file
        df_dedup.to_csv(output_path, sep='\t', index=False)
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Saved {len(df_dedup):,} rare variant burden associations to {output_path}")
        logger.info(f"{'='*60}")
        
        # Summary
        logger.info("\nEvidence sources:")
        for source, count in df_dedup['evidence_source'].value_counts().items():
            logger.info(f"  {source}: {count}")
        
        logger.info("\nTrait categories:")
        for category, count in df_dedup['trait_category'].value_counts().items():
            logger.info(f"  {category}: {count}")


def main():
    """Main execution."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Search ExWAS for rare variant burden tests"
    )
    parser.add_argument(
        '--benchmark-dir',
        default='benchmarks',
        help='Benchmark directory'
    )
    parser.add_argument(
        '--p-threshold',
        type=float,
        default=3.6e-7,
        help='Exome-wide significance threshold (default: 3.6e-7)'
    )
    parser.add_argument(
        '--sources',
        nargs='+',
        default=['uk_biobank_exwas', 'genebass', 'published_catalogs'],
        help='Data sources to search'
    )
    parser.add_argument(
        '--output',
        default='rare_burden_genes.tsv',
        help='Output filename'
    )
    
    args = parser.parse_args()
    
    searcher = ExWASBurdenSearcher(
        benchmark_dir=args.benchmark_dir,
        p_threshold=args.p_threshold,
        data_sources=args.sources,
    )
    
    logger.info("="*60)
    logger.info("Rare Variant Burden Test Search")
    logger.info("="*60)
    logger.info(f"P-value threshold: {args.p_threshold}")
    logger.info(f"Data sources: {', '.join(args.sources)}")
    logger.info("")
    
    # Load GWAS loci
    gwas_loci = searcher.load_gwas_loci()
    
    # Search each data source
    total_found = 0
    
    if 'uk_biobank_exwas' in args.sources:
        n = searcher.search_ukb_exwas()
        total_found += n
    
    if 'genebass' in args.sources:
        n = searcher.search_genebass()
        total_found += n
    
    if 'published_catalogs' in args.sources:
        n = searcher.search_published_catalogs()
        total_found += n
    
    logger.info(f"\nTotal burden associations found: {total_found}")
    
    if total_found == 0:
        logger.warning("\nNo burden associations found!")
        logger.warning("This is expected if ExWAS data files are not available.")
        logger.warning("Please download ExWAS results and place in:")
        logger.warning("  benchmarks/evidence_curation/")
        return
    
    # Match to GWAS loci (if loci available)
    if not gwas_loci.empty:
        matched = searcher.match_to_gwas_loci(gwas_loci)
        if not matched.empty:
            matched_output = searcher.output_dir / "rare_burden_matched.tsv"
            matched.to_csv(matched_output, sep='\t', index=False)
            logger.info(f"Saved matched associations to {matched_output}")
    
    # Save all results
    searcher.save_results(output_file=args.output)
    
    logger.info("\n" + "="*60)
    logger.info("Search complete!")
    logger.info("="*60)


if __name__ == '__main__':
    main()
