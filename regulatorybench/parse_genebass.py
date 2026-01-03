#!/usr/bin/env python3
"""
Extract rare variant burden genes from GeneBass UK Biobank exome-wide association study.

Criteria:
- p-value < 3.6e-7 (exome-wide significance threshold)
- Match phenotypes to GWAS trait ontology

Author: Generated for RegulatoryBench v3 platinum
Date: 2025-12-20
"""

import pandas as pd
import logging
from pathlib import Path
from typing import Dict, List, Set

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class GeneBassParser:
    """Parse GeneBass gene-level burden test results."""
    
    # GWAS trait ontology - same as mendelian parser
    TRAIT_ONTOLOGY = {
        'lipid': ['cholesterol', 'LDL', 'HDL', 'triglyceride', 'lipid', 'lipoprotein',
                  'hyperlipidemia', 'dyslipidemia'],
        'cardiovascular': ['coronary', 'myocardial', 'heart', 'cardiac', 'vascular',
                          'atherosclerosis', 'CAD', 'cardiomyopathy', 'arrhythmia',
                          'atrial fibrillation'],
        'diabetes': ['diabetes', 'T2D', 'glucose', 'insulin', 'glycemic'],
        'obesity': ['obesity', 'BMI', 'body mass', 'weight'],
        'blood_cells': ['erythrocyte', 'platelet', 'leukocyte', 'hemoglobin', 'RBC', 'WBC',
                       'blood count', 'hematocrit'],
        'immune': ['immune', 'autoimmune', 'inflammatory', 'rheumatoid', 'lupus', 'IBD',
                   'Crohn', 'ulcerative colitis', 'asthma'],
        'bone': ['bone', 'osteoporosis', 'skeletal', 'BMD', 'fracture'],
        'kidney': ['kidney', 'renal', 'nephropathy', 'CKD'],
        'liver': ['liver', 'hepatic', 'cirrhosis', 'NAFLD', 'hepatitis', 'ALT', 'AST'],
        'blood_pressure': ['hypertension', 'blood pressure', 'BP', 'systolic', 'diastolic'],
        'neurological': ['alzheimer', 'parkinson', 'dementia', 'epilepsy'],
    }
    
    # UK Biobank phenotype mappings
    UKBB_PHENOTYPE_MAP = {
        # Lipids
        '30780': 'lipid',  # LDL direct
        '30760': 'lipid',  # HDL cholesterol
        '30690': 'lipid',  # Total cholesterol
        '30870': 'lipid',  # Triglycerides
        
        # Cardiovascular
        'I25': 'cardiovascular',  # Chronic ischemic heart disease
        'I21': 'cardiovascular',  # Myocardial infarction
        'I48': 'cardiovascular',  # Atrial fibrillation
        'I50': 'cardiovascular',  # Heart failure
        
        # Diabetes
        'E11': 'diabetes',  # Type 2 diabetes
        '30740': 'diabetes',  # Glucose
        '30750': 'diabetes',  # HbA1c
        
        # Blood pressure
        '4080': 'blood_pressure',  # Systolic blood pressure
        '4079': 'blood_pressure',  # Diastolic blood pressure
        
        # Blood cells
        '30020': 'blood_cells',  # Red blood cell count
        '30080': 'blood_cells',  # Platelet count
        '30000': 'blood_cells',  # White blood cell count
        '30010': 'blood_cells',  # Lymphocyte count
        '30040': 'blood_cells',  # Neutrophil count
        
        # Kidney
        '30700': 'kidney',  # Creatinine
        'N18': 'kidney',  # Chronic kidney disease
        
        # Liver
        '30620': 'liver',  # ALT
        '30650': 'liver',  # AST
        'K74': 'liver',  # Fibrosis and cirrhosis
        
        # Obesity
        '21001': 'obesity',  # BMI
        
        # Immune
        'M05': 'immune',  # Rheumatoid arthritis
        'K50': 'immune',  # Crohn's disease
        'K51': 'immune',  # Ulcerative colitis
    }
    
    def __init__(self, data_dir: str = "data/external/genebass", output_dir: str = "benchmarks/evidence_curation"):
        self.data_dir = Path(data_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def match_phenotype(self, phenotype_str: str, phenotype_code: str = "") -> List[str]:
        """
        Match UK Biobank phenotype to GWAS trait categories.
        
        Parameters
        ----------
        phenotype_str : str
            Phenotype description
        phenotype_code : str
            UK Biobank phenotype code
            
        Returns
        -------
        list
            Matched GWAS trait categories
        """
        matched = []
        
        # Check UK Biobank code mapping
        if phenotype_code:
            for code, trait in self.UKBB_PHENOTYPE_MAP.items():
                if code in phenotype_code:
                    matched.append(trait)
        
        # Check phenotype text
        phenotype_lower = phenotype_str.lower()
        for trait_category, keywords in self.TRAIT_ONTOLOGY.items():
            if any(keyword.lower() in phenotype_lower for keyword in keywords):
                if trait_category not in matched:
                    matched.append(trait_category)
        
        return matched
    
    def parse_genebass(self, p_threshold: float = 3.6e-7) -> pd.DataFrame:
        """
        Parse GeneBass gene-level burden test results.
        
        Parameters
        ----------
        p_threshold : float
            P-value threshold (exome-wide significance, default 3.6e-7)
            
        Returns
        -------
        DataFrame
            Extracted rare burden genes with GWAS trait matches
        """
        logger.info("Parsing GeneBass gene-level burden test results...")
        
        genebass_file = self.data_dir / "gene_results.parquet"
        
        if not genebass_file.exists():
            logger.error(f"GeneBass file not found: {genebass_file}")
            logger.error("Please download from: https://genebass.org/downloads/gene_results.parquet")
            return pd.DataFrame()
        
        try:
            # Load GeneBass results
            logger.info("Loading GeneBass data...")
            df = pd.read_parquet(genebass_file)
            
            logger.info(f"Loaded {len(df):,} gene-phenotype associations")
            logger.info(f"Columns: {list(df.columns)}")
            
            # Filter for exome-wide significant results
            if 'Pvalue' in df.columns:
                df_sig = df[df['Pvalue'] < p_threshold].copy()
            elif 'pvalue' in df.columns:
                df_sig = df[df['pvalue'] < p_threshold].copy()
            elif 'P' in df.columns:
                df_sig = df[df['P'] < p_threshold].copy()
            else:
                logger.error(f"No p-value column found! Available columns: {list(df.columns)}")
                return pd.DataFrame()
            
            logger.info(f"Found {len(df_sig):,} exome-wide significant associations (p < {p_threshold:.2e})")
            
            if len(df_sig) == 0:
                logger.warning("No significant associations found!")
                return pd.DataFrame()
            
            # Extract gene symbols and phenotypes
            genes_curated = []
            
            for _, row in df_sig.iterrows():
                gene = str(row.get('gene_symbol', row.get('gene', row.get('GENE', '')))).strip()
                phenotype = str(row.get('phenotype', row.get('trait', row.get('TRAIT', '')))).strip()
                phenotype_code = str(row.get('phenocode', row.get('pheno_code', ''))).strip()
                pval = float(row.get('Pvalue', row.get('pvalue', row.get('P', 1.0))))
                
                if not gene or gene == 'nan':
                    continue
                
                # Match phenotype to GWAS traits
                matched_traits = self.match_phenotype(phenotype, phenotype_code)
                
                if matched_traits:
                    genes_curated.append({
                        'gene_symbol': gene,
                        'evidence_type': 'rare_variant_burden',
                        'evidence_source': 'GeneBass_UKB',
                        'phenotype': phenotype[:200],
                        'phenotype_code': phenotype_code if phenotype_code != 'nan' else '',
                        'matched_traits': '|'.join(matched_traits),
                        'p_value': f"{pval:.2e}",
                        'confidence': 'high',
                        'source_url': f"https://app.genebass.org/gene/{gene}",
                        'database_version': 'GeneBass_2023-11',
                        'date_curated': '2025-12-20'
                    })
            
            logger.info(f"Extracted {len(genes_curated)} genes matching GWAS traits")
            
            return pd.DataFrame(genes_curated)
            
        except Exception as e:
            logger.error(f"Error parsing GeneBass: {e}", exc_info=True)
            return pd.DataFrame()
    
    def deduplicate_and_save(self, df: pd.DataFrame, output_file: str = "rare_burden_genes.tsv"):
        """
        Deduplicate genes and save to file.
        
        Parameters
        ----------
        df : DataFrame
            Genes with burden test evidence
        output_file : str
            Output filename
        """
        if df.empty:
            logger.warning("No genes to save!")
            return
        
        logger.info("Deduplicating and aggregating evidence...")
        
        # Aggregate by gene
        gene_groups = df.groupby('gene_symbol')
        
        deduplicated = []
        for gene, group in gene_groups:
            # Combine phenotypes and p-values
            phenotypes = '; '.join(group['phenotype'].unique()[:5])  # Top 5
            best_pval = group['p_value'].min()
            
            # Combine matched traits
            all_traits = set()
            for traits_str in group['matched_traits']:
                all_traits.update(traits_str.split('|'))
            
            # Get representative row
            best_row = group.iloc[0].to_dict()
            best_row['phenotype'] = phenotypes
            best_row['matched_traits'] = '|'.join(sorted(all_traits))
            best_row['p_value'] = best_pval
            best_row['n_phenotypes'] = len(group)
            
            deduplicated.append(best_row)
        
        df_final = pd.DataFrame(deduplicated)
        
        # Sort by number of phenotypes and p-value
        df_final = df_final.sort_values(['n_phenotypes', 'p_value'], ascending=[False, True])
        
        # Save
        output_path = self.output_dir / output_file
        df_final.to_csv(output_path, sep='\t', index=False)
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Saved {len(df_final)} unique rare burden genes to {output_path}")
        logger.info(f"{'='*60}")
        
        # Summary statistics
        logger.info("\nSummary by GWAS trait category:")
        all_trait_counts = {}
        for traits_str in df_final['matched_traits']:
            for trait in traits_str.split('|'):
                all_trait_counts[trait] = all_trait_counts.get(trait, 0) + 1
        
        for trait, count in sorted(all_trait_counts.items(), key=lambda x: -x[1])[:15]:
            logger.info(f"  {trait}: {count}")
        
        logger.info(f"\nTop genes by number of phenotypes:")
        for _, row in df_final.head(10).iterrows():
            logger.info(f"  {row['gene_symbol']}: {row['n_phenotypes']} phenotypes ({row['matched_traits']})")


def main():
    """Main execution."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Parse GeneBass rare variant burden test results"
    )
    parser.add_argument(
        '--data-dir',
        default='data/external/genebass',
        help='GeneBass data directory'
    )
    parser.add_argument(
        '--output-dir',
        default='benchmarks/evidence_curation',
        help='Output directory'
    )
    parser.add_argument(
        '--p-threshold',
        type=float,
        default=3.6e-7,
        help='P-value threshold (exome-wide significance, default: 3.6e-7)'
    )
    parser.add_argument(
        '--output',
        default='rare_burden_genes.tsv',
        help='Output filename'
    )
    
    args = parser.parse_args()
    
    parser_obj = GeneBassParser(
        data_dir=args.data_dir,
        output_dir=args.output_dir
    )
    
    logger.info("="*60)
    logger.info("GeneBass Rare Variant Burden Gene Curation")
    logger.info("="*60)
    logger.info("")
    
    # Parse GeneBass
    df = parser_obj.parse_genebass(p_threshold=args.p_threshold)
    
    if df.empty:
        logger.error("\nNo genes extracted!")
        logger.error("Please ensure GeneBass data is downloaded:")
        logger.error("  - data/external/genebass/gene_results.parquet")
        return 1
    
    # Deduplicate and save
    parser_obj.deduplicate_and_save(df, output_file=args.output)
    
    logger.info("\n" + "="*60)
    logger.info("Curation complete!")
    logger.info("="*60)
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
