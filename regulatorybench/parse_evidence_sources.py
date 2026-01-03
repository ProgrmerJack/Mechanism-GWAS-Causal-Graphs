#!/usr/bin/env python3
"""
Parse External Evidence Sources for v3 Platinum Gold Standard

Parses downloaded ClinVar and Gene2Phenotype data to extract Mendelian disease genes.

Author: Generated for RegulatoryBench v3 platinum
Date: 2025-12-20
"""

import pandas as pd
import gzip
import logging
from pathlib import Path
from typing import Dict, List, Set
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class EvidenceSourceParser:
    """Parse external evidence sources for Mendelian disease genes."""
    
    # GWAS trait ontology for matching Mendelian diseases to GWAS traits
    TRAIT_ONTOLOGY = {
        'lipid': ['cholesterol', 'LDL', 'HDL', 'triglyceride', 'lipid', 'lipoprotein',
                  'hyperlipidemia', 'dyslipidemia', 'familial hypercholesterolemia',
                  'hypercholesterolemia', 'sitosterolemia', 'abetalipoproteinemia'],
        'cardiovascular': ['coronary', 'myocardial', 'heart', 'cardiac', 'vascular',
                          'atherosclerosis', 'CAD', 'cardiomyopathy', 'arrhythmia',
                          'atrial fibrillation', 'QT', 'Brugada'],
        'diabetes': ['diabetes', 'T2D', 'glucose', 'insulin', 'glycemic', 'MODY',
                     'hyperinsulinism', 'hypoglycemia', 'maturity-onset diabetes'],
        'obesity': ['obesity', 'BMI', 'body mass', 'weight', 'adipose', 'leptin'],
        'blood_cells': ['erythrocyte', 'platelet', 'leukocyte', 'hemoglobin', 'RBC', 'WBC',
                       'anemia', 'thalassemia', 'sickle cell', 'hematologic',
                       'thrombocytopenia', 'neutropenia'],
        'immune': ['immune', 'autoimmune', 'inflammatory', 'rheumatoid', 'lupus', 'IBD',
                   'Crohn', 'ulcerative colitis', 'immunodeficiency', 'celiac'],
        'bone': ['bone', 'osteoporosis', 'skeletal', 'BMD', 'fracture', 'osteogenesis',
                 'osteopetrosis', 'rickets'],
        'kidney': ['kidney', 'renal', 'nephropathy', 'CKD', 'nephritis', 'nephrotic'],
        'liver': ['liver', 'hepatic', 'cirrhosis', 'NAFLD', 'hepatitis', 'cholestasis'],
        'blood_pressure': ['hypertension', 'blood pressure', 'BP', 'hypotension'],
        'neurological': ['alzheimer', 'parkinson', 'dementia', 'epilepsy', 'ataxia',
                        'neuropathy', 'leukodystrophy'],
    }
    
    def __init__(self, data_dir: str = "data/external", output_dir: str = "benchmarks/evidence_curation"):
        self.data_dir = Path(data_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.mendelian_genes = []
        
    def match_phenotype_to_traits(self, phenotype: str) -> List[str]:
        """
        Match phenotype string to GWAS trait categories.
        
        Parameters
        ----------
        phenotype : str
            Disease/phenotype description
            
        Returns
        -------
        list
            Matched GWAS trait categories
        """
        phenotype_lower = phenotype.lower()
        matched = []
        
        for trait_category, keywords in self.TRAIT_ONTOLOGY.items():
            if any(keyword.lower() in phenotype_lower for keyword in keywords):
                matched.append(trait_category)
        
        return matched
    
    def parse_clinvar(self, min_pathogenic: int = 2) -> int:
        """
        Parse ClinVar variant_summary.txt.gz for pathogenic variants.
        
        Parameters
        ----------
        min_pathogenic : int
            Minimum pathogenic variants per gene
            
        Returns
        -------
        int
            Number of genes extracted
        """
        logger.info("Parsing ClinVar variant summary...")
        
        clinvar_file = self.data_dir / "clinvar" / "variant_summary.txt.gz"
        
        if not clinvar_file.exists():
            logger.error(f"ClinVar file not found: {clinvar_file}")
            logger.error("Please download from: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz")
            return 0
        
        try:
            # Read ClinVar data (this is a large file, ~3GB)
            logger.info("Reading ClinVar data (this may take a few minutes)...")
            
            with gzip.open(clinvar_file, 'rt') as f:
                # Read header
                header = f.readline().strip().split('\t')
                logger.info(f"ClinVar columns: {len(header)}")
                
                # Find column indices
                gene_idx = header.index('GeneSymbol')
                sig_idx = header.index('ClinicalSignificance')
                pheno_idx = header.index('PhenotypeList')
                assembly_idx = header.index('Assembly')
                review_idx = header.index('ReviewStatus')
                
                # Parse variants
                gene_data = defaultdict(lambda: {'pathogenic_count': 0, 'phenotypes': set()})
                
                for i, line in enumerate(f):
                    if (i + 1) % 100000 == 0:
                        logger.info(f"  Processed {i+1:,} variants...")
                    
                    fields = line.strip().split('\t')
                    
                    # Skip malformed lines (incomplete downloads)
                    if len(fields) < max(gene_idx, sig_idx, pheno_idx, assembly_idx, review_idx) + 1:
                        continue
                    
                    # Check assembly
                    if fields[assembly_idx] != 'GRCh38':
                        continue
                    
                    gene = fields[gene_idx]
                    sig = fields[sig_idx].lower()
                    pheno = fields[pheno_idx]
                    
                    # Skip if no gene or conflicting/uncertain
                    if not gene or 'conflicting' in sig or 'uncertain' in sig:
                        continue
                    
                    # Check for pathogenic/likely pathogenic
                    if 'pathogenic' in sig or 'likely pathogenic' in sig:
                        gene_data[gene]['pathogenic_count'] += 1
                        if pheno and pheno != '-':
                            gene_data[gene]['phenotypes'].update(pheno.split(';'))
            
            logger.info(f"Parsed ClinVar: {len(gene_data):,} genes with pathogenic variants")
            
            # Filter and match to GWAS traits
            count = 0
            for gene, data in gene_data.items():
                if data['pathogenic_count'] < min_pathogenic:
                    continue
                
                # Match phenotypes to GWAS traits
                phenotype_str = ' '.join(data['phenotypes'])
                matched_traits = self.match_phenotype_to_traits(phenotype_str)
                
                if matched_traits:
                    self.mendelian_genes.append({
                        'gene_symbol': gene,
                        'evidence_type': 'mendelian_disease',
                        'evidence_source': 'ClinVar_pathogenic',
                        'n_pathogenic_variants': data['pathogenic_count'],
                        'matched_traits': '|'.join(matched_traits),
                        'phenotypes': phenotype_str[:500] if phenotype_str else '',
                        'confidence': 'very_high' if data['pathogenic_count'] >= 10 else 'high',
                        'source_url': f"https://www.ncbi.nlm.nih.gov/clinvar/?term={gene}[gene]",
                        'database_version': 'ClinVar_2025-12',
                        'date_curated': '2025-12-20'
                    })
                    count += 1
            
            logger.info(f"Extracted {count} genes from ClinVar matching GWAS traits")
            return count
            
        except Exception as e:
            logger.error(f"Error parsing ClinVar: {e}", exc_info=True)
            return 0
    
    def parse_gene2phenotype(self) -> int:
        """
        Parse Gene2Phenotype panel files.
        
        Returns
        -------
        int
            Number of genes extracted
        """
        logger.info("Parsing Gene2Phenotype panels...")
        
        g2p_dir = self.data_dir / "g2p"
        
        # Find all downloaded panel files
        panel_files = list(g2p_dir.glob("*.csv"))
        
        if not panel_files:
            logger.warning(f"No Gene2Phenotype files found in {g2p_dir}")
            return 0
        
        logger.info(f"Found {len(panel_files)} Gene2Phenotype panel files")
        
        count = 0
        for panel_file in panel_files:
            try:
                logger.info(f"  Parsing {panel_file.name}...")
                
                # Read panel file
                df = pd.read_csv(panel_file)
                
                logger.info(f"    Loaded {len(df)} entries")
                
                # Extract genes with confident associations
                for _, row in df.iterrows():
                    gene = str(row.get('gene symbol', row.get('gene_symbol', ''))).strip()
                    disease = str(row.get('disease name', row.get('disease_name', ''))).strip()
                    confidence = str(row.get('confidence category', row.get('confidence', 'unknown'))).lower()
                    
                    if not gene or not disease or gene == 'nan' or disease == 'nan':
                        continue
                    
                    # Filter for confident associations
                    if 'confirmed' not in confidence and 'definitive' not in confidence:
                        continue
                    
                    # Match disease to GWAS traits
                    matched_traits = self.match_phenotype_to_traits(disease)
                    
                    if matched_traits:
                        self.mendelian_genes.append({
                            'gene_symbol': gene,
                            'evidence_type': 'mendelian_disease',
                            'evidence_source': f"Gene2Phenotype_{panel_file.stem.replace('.csv', '')}",
                            'disease_name': disease[:200],
                            'matched_traits': '|'.join(matched_traits),
                            'confidence': 'very_high',
                            'source_url': 'https://www.ebi.ac.uk/gene2phenotype',
                            'database_version': 'G2P_2025-12',
                            'date_curated': '2025-12-20'
                        })
                        count += 1
                
            except Exception as e:
                logger.error(f"Error parsing {panel_file}: {e}")
                continue
        
        logger.info(f"Extracted {count} genes from Gene2Phenotype matching GWAS traits")
        return count
    
    def deduplicate_and_save(self, output_file: str = "mendelian_genes_curated.tsv"):
        """
        Deduplicate genes and save to file.
        
        Parameters
        ----------
        output_file : str
            Output filename
        """
        if not self.mendelian_genes:
            logger.warning("No Mendelian genes to save!")
            return
        
        logger.info("Deduplicating and aggregating evidence...")
        
        df = pd.DataFrame(self.mendelian_genes)
        
        # Aggregate by gene
        gene_groups = df.groupby('gene_symbol')
        
        deduplicated = []
        for gene, group in gene_groups:
            # Combine evidence sources
            sources = '|'.join(sorted(group['evidence_source'].unique()))
            
            # Combine matched traits
            all_traits = set()
            for traits_str in group['matched_traits']:
                all_traits.update(traits_str.split('|'))
            
            # Highest confidence
            confidences = group['confidence'].tolist()
            best_confidence = 'very_high' if 'very_high' in confidences else 'high'
            
            # Get representative row
            best_row = group.iloc[0].to_dict()
            best_row['evidence_source'] = sources
            best_row['matched_traits'] = '|'.join(sorted(all_traits))
            best_row['confidence'] = best_confidence
            best_row['n_sources'] = len(group)
            
            deduplicated.append(best_row)
        
        df_final = pd.DataFrame(deduplicated)
        
        # Sort by confidence and number of sources
        df_final['conf_rank'] = df_final['confidence'].map({'very_high': 1, 'high': 2, 'medium': 3})
        df_final = df_final.sort_values(['conf_rank', 'n_sources'], ascending=[True, False])
        df_final = df_final.drop('conf_rank', axis=1)
        
        # Save
        output_path = self.output_dir / output_file
        df_final.to_csv(output_path, sep='\t', index=False)
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Saved {len(df_final)} unique Mendelian genes to {output_path}")
        logger.info(f"{'='*60}")
        
        # Summary statistics
        logger.info("\nSummary by evidence source:")
        for source, count in df['evidence_source'].value_counts().items():
            logger.info(f"  {source}: {count}")
        
        logger.info("\nSummary by GWAS trait category:")
        all_trait_counts = defaultdict(int)
        for traits_str in df_final['matched_traits']:
            for trait in traits_str.split('|'):
                all_trait_counts[trait] += 1
        
        for trait, count in sorted(all_trait_counts.items(), key=lambda x: -x[1])[:15]:
            logger.info(f"  {trait}: {count}")
        
        logger.info(f"\nConfidence breakdown:")
        for conf, count in df_final['confidence'].value_counts().items():
            logger.info(f"  {conf}: {count}")


def main():
    """Main execution."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Parse external evidence sources for Mendelian disease genes"
    )
    parser.add_argument(
        '--data-dir',
        default='data/external',
        help='External data directory'
    )
    parser.add_argument(
        '--output-dir',
        default='benchmarks/evidence_curation',
        help='Output directory'
    )
    parser.add_argument(
        '--min-pathogenic',
        type=int,
        default=2,
        help='Minimum pathogenic variants per gene in ClinVar (default: 2)'
    )
    parser.add_argument(
        '--output',
        default='mendelian_genes_curated.tsv',
        help='Output filename'
    )
    
    args = parser.parse_args()
    
    parser_obj = EvidenceSourceParser(
        data_dir=args.data_dir,
        output_dir=args.output_dir
    )
    
    logger.info("="*60)
    logger.info("Mendelian Disease Gene Curation")
    logger.info("="*60)
    logger.info("")
    
    # Parse ClinVar
    n_clinvar = parser_obj.parse_clinvar(min_pathogenic=args.min_pathogenic)
    
    # Parse Gene2Phenotype
    n_g2p = parser_obj.parse_gene2phenotype()
    
    total = n_clinvar + n_g2p
    
    if total == 0:
        logger.error("\nNo genes extracted from any source!")
        logger.error("Please ensure data files are downloaded:")
        logger.error("  - data/external/clinvar/variant_summary.txt.gz")
        logger.error("  - data/external/g2p/*.csv.gz")
        return 1
    
    # Deduplicate and save
    parser_obj.deduplicate_and_save(output_file=args.output)
    
    logger.info("\n" + "="*60)
    logger.info("Curation complete!")
    logger.info("="*60)
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
