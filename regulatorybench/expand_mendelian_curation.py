#!/usr/bin/env python3
"""
Expand Mendelian Disease Gene Curation

Comprehensive curation of Mendelian disease genes from multiple sources:
- ClinVar (pathogenic/likely pathogenic variants)
- OMIM (disease genes)
- Gene2Phenotype (curated panels)
- ClinGen (gene-disease validity)
- OpenTargets Genetics (Mendelian subset)

Goal: Expand from current 44 genes to 100-150 genes with full provenance.

Author: Generated for RegulatoryBench v3 platinum
Date: 2025-01
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Set
import requests
import json
from datetime import datetime

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class MendelianCurator:
    """Curate Mendelian disease genes from multiple authoritative sources."""
    
    def __init__(self, benchmark_dir: str = "benchmarks"):
        self.benchmark_dir = Path(benchmark_dir)
        self.output_dir = self.benchmark_dir / "evidence_curation"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load GWAS loci for trait matching
        self.gwas_loci = None
        self.trait_ontology = self._load_trait_ontology()
        
        # Curated genes storage
        self.curated_genes = []
        
    def _load_trait_ontology(self) -> Dict[str, List[str]]:
        """Load GWAS trait ontology for matching Mendelian diseases."""
        ontology = {
            # Lipid/cardiovascular
            'lipid': ['cholesterol', 'LDL', 'HDL', 'triglyceride', 'lipid', 'lipoprotein', 
                      'hyperlipidemia', 'dyslipidemia', 'familial hypercholesterolemia'],
            'cardiovascular': ['coronary', 'myocardial', 'heart', 'cardiac', 'vascular', 
                               'atherosclerosis', 'CAD'],
            
            # Metabolic
            'diabetes': ['diabetes', 'T2D', 'glucose', 'insulin', 'glycemic', 'MODY'],
            'obesity': ['obesity', 'BMI', 'body mass', 'weight', 'adipose'],
            
            # Blood/immune
            'blood_cells': ['erythrocyte', 'platelet', 'leukocyte', 'hemoglobin', 'RBC', 'WBC',
                            'anemia', 'thalassemia', 'sickle cell'],
            'immune': ['immune', 'autoimmune', 'inflammatory', 'rheumatoid', 'lupus', 'IBD',
                       'Crohn', 'ulcerative colitis'],
            
            # Other
            'bone': ['bone', 'osteoporosis', 'skeletal', 'BMD', 'fracture'],
            'kidney': ['kidney', 'renal', 'nephropathy', 'CKD'],
            'liver': ['liver', 'hepatic', 'cirrhosis', 'NAFLD'],
        }
        return ontology
    
    def curate_from_clinvar(self, min_pathogenic: int = 3) -> int:
        """
        Curate genes from ClinVar with pathogenic/likely pathogenic variants.
        
        Uses ClinVar variant summary to identify genes with multiple
        high-confidence pathogenic variants in traits matching GWAS.
        
        Parameters
        ----------
        min_pathogenic : int
            Minimum number of pathogenic variants required per gene.
            
        Returns
        -------
        int
            Number of genes curated.
        """
        logger.info("Curating genes from ClinVar...")
        
        # Download ClinVar variant summary (note: ~3GB file)
        # https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
        clinvar_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
        
        logger.info(f"Downloading ClinVar variant summary from {clinvar_url}")
        logger.warning("This may take several minutes (file is ~3GB)...")
        
        try:
            # Download and read ClinVar data
            clinvar_df = pd.read_csv(
                clinvar_url,
                sep='\t',
                compression='gzip',
                usecols=['GeneSymbol', 'ClinicalSignificance', 'PhenotypeList', 
                         'ReviewStatus', 'Assembly', 'Chromosome', 'VariationID'],
                low_memory=False
            )
            
            logger.info(f"Loaded {len(clinvar_df):,} ClinVar variants")
            
            # Filter for GRCh38, pathogenic/likely pathogenic
            clinvar_df = clinvar_df[
                (clinvar_df['Assembly'] == 'GRCh38') &
                (clinvar_df['ClinicalSignificance'].str.contains(
                    'Pathogenic|Likely pathogenic', case=False, na=False
                )) &
                (~clinvar_df['ClinicalSignificance'].str.contains(
                    'Conflicting|Uncertain', case=False, na=False
                ))
            ]
            
            logger.info(f"Filtered to {len(clinvar_df):,} pathogenic/likely pathogenic variants")
            
            # Count pathogenic variants per gene
            gene_counts = clinvar_df.groupby('GeneSymbol').size().reset_index(name='n_pathogenic')
            gene_counts = gene_counts[gene_counts['n_pathogenic'] >= min_pathogenic]
            
            logger.info(f"Found {len(gene_counts):,} genes with >={min_pathogenic} pathogenic variants")
            
            # Match phenotypes to GWAS traits
            count = 0
            for _, row in gene_counts.iterrows():
                gene = row['GeneSymbol']
                
                # Get phenotypes for this gene
                gene_variants = clinvar_df[clinvar_df['GeneSymbol'] == gene]
                phenotypes = ' '.join(gene_variants['PhenotypeList'].fillna('').values)
                
                # Match to GWAS trait ontology
                matched_traits = self._match_phenotype_to_gwas(phenotypes)
                
                if matched_traits:
                    self.curated_genes.append({
                        'gene_symbol': gene,
                        'evidence_type': 'mendelian_disease',
                        'evidence_source': 'ClinVar',
                        'source_url': f'https://www.ncbi.nlm.nih.gov/clinvar/?term={gene}[gene]',
                        'n_pathogenic_variants': int(row['n_pathogenic']),
                        'matched_traits': ', '.join(matched_traits),
                        'phenotypes': phenotypes[:200],  # Truncate
                        'confidence': 'very_high' if row['n_pathogenic'] >= 10 else 'high',
                        'date_curated': datetime.now().strftime('%Y-%m-%d')
                    })
                    count += 1
            
            logger.info(f"Curated {count} genes from ClinVar matching GWAS traits")
            return count
            
        except Exception as e:
            logger.error(f"Error curating from ClinVar: {e}")
            logger.warning("Continuing with other sources...")
            return 0
    
    def curate_from_omim(self, omim_api_key: str = None) -> int:
        """
        Curate genes from OMIM disease associations.
        
        Requires OMIM API key (free registration at https://omim.org/api).
        Alternatively, can use genemap2.txt if downloaded manually.
        
        Parameters
        ----------
        omim_api_key : str, optional
            OMIM API key for programmatic access.
            
        Returns
        -------
        int
            Number of genes curated.
        """
        logger.info("Curating genes from OMIM...")
        
        # Check for manually downloaded genemap2.txt
        genemap_file = Path('data/external/genemap2.txt')
        
        if genemap_file.exists():
            logger.info(f"Using OMIM genemap2.txt from {genemap_file}")
            return self._curate_from_genemap(genemap_file)
        
        elif omim_api_key:
            logger.info("Using OMIM API...")
            return self._curate_from_omim_api(omim_api_key)
        
        else:
            logger.warning("""
OMIM data not available. To curate from OMIM:

Option 1: Manual download
  1. Register at https://omim.org/downloads
  2. Download genemap2.txt
  3. Place in data/external/genemap2.txt

Option 2: API access
  1. Register at https://omim.org/api
  2. Get API key
  3. Run: expand_mendelian_curation.py --omim-api-key YOUR_KEY

Continuing with other sources...
""")
            return 0
    
    def _curate_from_genemap(self, genemap_file: Path) -> int:
        """Curate from OMIM genemap2.txt file."""
        try:
            # Read genemap2.txt (tab-delimited)
            df = pd.read_csv(genemap_file, sep='\t', comment='#', low_memory=False)
            
            logger.info(f"Loaded {len(df):,} OMIM gene-disease associations")
            
            count = 0
            for _, row in df.iterrows():
                gene = str(row.get('Gene Symbols', '')).split(',')[0].strip()
                phenotypes = str(row.get('Phenotypes', ''))
                mim_number = str(row.get('MIM Number', ''))
                
                if not gene or not phenotypes or phenotypes == 'nan':
                    continue
                
                # Match phenotypes to GWAS traits
                matched_traits = self._match_phenotype_to_gwas(phenotypes)
                
                if matched_traits:
                    self.curated_genes.append({
                        'gene_symbol': gene,
                        'evidence_type': 'mendelian_disease',
                        'evidence_source': 'OMIM',
                        'source_url': f'https://omim.org/entry/{mim_number}',
                        'omim_id': mim_number,
                        'matched_traits': ', '.join(matched_traits),
                        'phenotypes': phenotypes[:200],
                        'confidence': 'very_high',
                        'date_curated': datetime.now().strftime('%Y-%m-%d')
                    })
                    count += 1
            
            logger.info(f"Curated {count} genes from OMIM matching GWAS traits")
            return count
            
        except Exception as e:
            logger.error(f"Error curating from OMIM genemap: {e}")
            return 0
    
    def _curate_from_omim_api(self, api_key: str) -> int:
        """Curate from OMIM API."""
        # Implementation for OMIM API would go here
        # Requires API key and handling of rate limits
        logger.warning("OMIM API curation not yet implemented. Use genemap2.txt instead.")
        return 0
    
    def curate_from_gene2phenotype(self) -> int:
        """
        Curate genes from Gene2Phenotype panels.
        
        Gene2Phenotype provides curated gene-disease associations for:
        - Developmental disorders (DDG2P)
        - Cancer
        - Cardiac disorders
        - Skin disorders
        - Eye disorders
        
        Returns
        -------
        int
            Number of genes curated.
        """
        logger.info("Curating genes from Gene2Phenotype...")
        
        # Gene2Phenotype downloads (CSV format)
        g2p_urls = {
            'DDG2P': 'https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz',
            'Cancer': 'https://www.ebi.ac.uk/gene2phenotype/downloads/CancerG2P.csv.gz',
            'Cardiac': 'https://www.ebi.ac.uk/gene2phenotype/downloads/CardiacG2P.csv.gz',
            'Skin': 'https://www.ebi.ac.uk/gene2phenotype/downloads/SkinG2P.csv.gz',
            'Eye': 'https://www.ebi.ac.uk/gene2phenotype/downloads/EyeG2P.csv.gz',
        }
        
        count = 0
        for panel_name, url in g2p_urls.items():
            try:
                logger.info(f"Downloading {panel_name} panel from {url}")
                df = pd.read_csv(url, compression='gzip')
                
                logger.info(f"Loaded {len(df):,} associations from {panel_name}")
                
                for _, row in df.iterrows():
                    gene = str(row.get('gene symbol', ''))
                    disease = str(row.get('disease name', ''))
                    confidence = str(row.get('confidence category', '')).lower()
                    
                    if not gene or not disease:
                        continue
                    
                    # Match disease to GWAS traits
                    matched_traits = self._match_phenotype_to_gwas(disease)
                    
                    if matched_traits:
                        # Map G2P confidence to our scale
                        conf_mapping = {
                            'definitive': 'very_high',
                            'strong': 'very_high',
                            'moderate': 'high',
                            'limited': 'medium',
                            'both dd and if': 'high',
                            'confirmed': 'very_high',
                            'probable': 'high',
                            'possible': 'medium',
                        }
                        
                        our_confidence = conf_mapping.get(confidence, 'medium')
                        
                        self.curated_genes.append({
                            'gene_symbol': gene,
                            'evidence_type': 'mendelian_disease',
                            'evidence_source': f'Gene2Phenotype_{panel_name}',
                            'source_url': f'https://www.ebi.ac.uk/gene2phenotype/search?panel={panel_name}&search_term={gene}',
                            'matched_traits': ', '.join(matched_traits),
                            'disease': disease[:200],
                            'g2p_confidence': confidence,
                            'confidence': our_confidence,
                            'date_curated': datetime.now().strftime('%Y-%m-%d')
                        })
                        count += 1
                
            except Exception as e:
                logger.warning(f"Could not download {panel_name} panel: {e}")
                continue
        
        logger.info(f"Curated {count} genes from Gene2Phenotype matching GWAS traits")
        return count
    
    def curate_from_clingen(self) -> int:
        """
        Curate genes from ClinGen gene-disease validity curations.
        
        ClinGen provides expert-curated gene-disease associations with
        evidence-based validity classifications.
        
        Returns
        -------
        int
            Number of genes curated.
        """
        logger.info("Curating genes from ClinGen...")
        
        # ClinGen gene-disease validity download
        clingen_url = "https://search.clinicalgenome.org/kb/downloads/gene-disease-validity-download.csv"
        
        try:
            logger.info(f"Downloading ClinGen curations from {clingen_url}")
            df = pd.read_csv(clingen_url)
            
            logger.info(f"Loaded {len(df):,} ClinGen gene-disease curations")
            
            # Filter for high-confidence curations
            df = df[df['CLASSIFICATION'].isin([
                'Definitive', 'Strong', 'Moderate', 'Definitive', 'Definitive '
            ])]
            
            logger.info(f"Filtered to {len(df):,} high-confidence curations")
            
            count = 0
            for _, row in df.iterrows():
                gene = str(row.get('GENE SYMBOL', ''))
                disease = str(row.get('DISEASE LABEL', ''))
                classification = str(row.get('CLASSIFICATION', ''))
                uuid = str(row.get('CLINGEN UUID', ''))
                
                if not gene or not disease:
                    continue
                
                # Match disease to GWAS traits
                matched_traits = self._match_phenotype_to_gwas(disease)
                
                if matched_traits:
                    # Map ClinGen classification to confidence
                    conf_mapping = {
                        'Definitive': 'very_high',
                        'Strong': 'very_high',
                        'Moderate': 'high',
                    }
                    
                    our_confidence = conf_mapping.get(classification, 'medium')
                    
                    self.curated_genes.append({
                        'gene_symbol': gene,
                        'evidence_type': 'mendelian_disease',
                        'evidence_source': 'ClinGen',
                        'source_url': f'https://search.clinicalgenome.org/kb/gene-validity/{uuid}',
                        'clingen_uuid': uuid,
                        'matched_traits': ', '.join(matched_traits),
                        'disease': disease[:200],
                        'clingen_classification': classification,
                        'confidence': our_confidence,
                        'date_curated': datetime.now().strftime('%Y-%m-%d')
                    })
                    count += 1
            
            logger.info(f"Curated {count} genes from ClinGen matching GWAS traits")
            return count
            
        except Exception as e:
            logger.error(f"Error curating from ClinGen: {e}")
            return 0
    
    def curate_from_opentargets(self) -> int:
        """
        Curate Mendelian genes from OpenTargets Genetics.
        
        Note: OpenTargets L2G scores are computational, so we ONLY use
        their manually curated Mendelian disease genes subset.
        
        Returns
        -------
        int
            Number of genes curated.
        """
        logger.info("Curating Mendelian genes from OpenTargets Genetics...")
        
        # OpenTargets GraphQL API for gene-disease associations
        # We query for genes with "genetic association" evidence from Mendelian sources
        
        opentargets_api = "https://api.platform.opentargets.org/api/v4/graphql"
        
        # This is a simplified example - full implementation would paginate through results
        query = """
        query mendelianGenes {
          diseases {
            rows {
              id
              name
              therapeuticAreas {
                id
                name
              }
            }
          }
        }
        """
        
        # Note: Full implementation would query for:
        # - Evidence with datatype "genetic_association"
        # - Evidence sources: "gene2phenotype", "orphanet", "genomics_england"
        # - Filter for high-confidence associations
        
        logger.warning("""
OpenTargets Genetics curation requires GraphQL API queries.
Full implementation would:
1. Query for genes with genetic_association evidence
2. Filter for Mendelian sources (Gene2Phenotype, Orphanet, ClinVar)
3. Match diseases to GWAS trait ontology
4. Exclude computational L2G scores (not independent evidence)

Skipping OpenTargets for now - using other authoritative sources instead.
""")
        
        return 0
    
    def _match_phenotype_to_gwas(self, phenotype_text: str) -> List[str]:
        """
        Match disease/phenotype text to GWAS trait categories.
        
        Parameters
        ----------
        phenotype_text : str
            Disease or phenotype description.
            
        Returns
        -------
        list
            List of matched GWAS trait categories.
        """
        if not isinstance(phenotype_text, str):
            return []
        
        phenotype_lower = phenotype_text.lower()
        matched = []
        
        for trait_category, keywords in self.trait_ontology.items():
            if any(keyword.lower() in phenotype_lower for keyword in keywords):
                matched.append(trait_category)
        
        return list(set(matched))  # Unique
    
    def deduplicate_and_prioritize(self) -> pd.DataFrame:
        """
        Deduplicate genes and prioritize by evidence strength.
        
        When a gene appears from multiple sources:
        1. Prioritize higher confidence
        2. Combine evidence sources
        3. Keep most informative metadata
        
        Returns
        -------
        pd.DataFrame
            Deduplicated curated genes.
        """
        logger.info("Deduplicating and prioritizing curated genes...")
        
        if not self.curated_genes:
            logger.warning("No genes curated!")
            return pd.DataFrame()
        
        df = pd.DataFrame(self.curated_genes)
        
        logger.info(f"Before deduplication: {len(df):,} gene-trait associations")
        
        # Group by gene and matched traits
        grouped = df.groupby(['gene_symbol', 'matched_traits']).agg({
            'evidence_source': lambda x: ', '.join(sorted(set(x))),
            'confidence': lambda x: self._highest_confidence(x),
            'date_curated': 'first',
            'source_url': 'first',  # Keep first URL
        }).reset_index()
        
        logger.info(f"After deduplication: {len(grouped):,} unique gene-trait pairs")
        
        # Count genes per trait category
        trait_counts = {}
        for traits in grouped['matched_traits']:
            for trait in traits.split(', '):
                trait_counts[trait] = trait_counts.get(trait, 0) + 1
        
        logger.info("\nGenes per trait category:")
        for trait, count in sorted(trait_counts.items(), key=lambda x: -x[1]):
            logger.info(f"  {trait}: {count}")
        
        return grouped
    
    def _highest_confidence(self, confidences: pd.Series) -> str:
        """Select highest confidence level from multiple values."""
        conf_order = ['very_high', 'high', 'medium', 'low']
        for conf in conf_order:
            if conf in confidences.values:
                return conf
        return 'medium'
    
    def save_results(self, output_file: str = "mendelian_genes_expanded.tsv"):
        """
        Save curated Mendelian genes to file.
        
        Parameters
        ----------
        output_file : str
            Output filename.
        """
        df = self.deduplicate_and_prioritize()
        
        if df.empty:
            logger.error("No genes to save!")
            return
        
        output_path = self.output_dir / output_file
        df.to_csv(output_path, sep='\t', index=False)
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Saved {len(df):,} curated Mendelian genes to {output_path}")
        logger.info(f"{'='*60}")
        
        # Summary statistics
        logger.info("\nSummary Statistics:")
        logger.info(f"  Total unique genes: {df['gene_symbol'].nunique()}")
        logger.info(f"  Very high confidence: {len(df[df['confidence']=='very_high'])}")
        logger.info(f"  High confidence: {len(df[df['confidence']=='high'])}")
        logger.info(f"  Medium confidence: {len(df[df['confidence']=='medium'])}")
        
        logger.info("\nEvidence sources used:")
        sources = df['evidence_source'].str.split(', ').explode().value_counts()
        for source, count in sources.items():
            logger.info(f"  {source}: {count}")


def main():
    """Main execution."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Expand Mendelian disease gene curation for RegulatoryBench v3"
    )
    parser.add_argument(
        '--benchmark-dir',
        default='benchmarks',
        help='Benchmark directory (default: benchmarks)'
    )
    parser.add_argument(
        '--omim-api-key',
        help='OMIM API key (optional, for programmatic access)'
    )
    parser.add_argument(
        '--output',
        default='mendelian_genes_expanded.tsv',
        help='Output filename (default: mendelian_genes_expanded.tsv)'
    )
    parser.add_argument(
        '--sources',
        nargs='+',
        choices=['clinvar', 'omim', 'gene2phenotype', 'clingen', 'all'],
        default=['all'],
        help='Evidence sources to use (default: all)'
    )
    
    args = parser.parse_args()
    
    curator = MendelianCurator(benchmark_dir=args.benchmark_dir)
    
    # Determine which sources to use
    sources = args.sources
    if 'all' in sources:
        sources = ['clinvar', 'omim', 'gene2phenotype', 'clingen']
    
    logger.info("="*60)
    logger.info("Mendelian Disease Gene Curation")
    logger.info("="*60)
    logger.info(f"Sources to use: {', '.join(sources)}")
    logger.info("")
    
    # Curate from each source
    total_curated = 0
    
    if 'clinvar' in sources:
        total_curated += curator.curate_from_clinvar(min_pathogenic=3)
    
    if 'omim' in sources:
        total_curated += curator.curate_from_omim(omim_api_key=args.omim_api_key)
    
    if 'gene2phenotype' in sources:
        total_curated += curator.curate_from_gene2phenotype()
    
    if 'clingen' in sources:
        total_curated += curator.curate_from_clingen()
    
    logger.info(f"\nTotal associations curated: {total_curated}")
    
    # Save results
    curator.save_results(output_file=args.output)
    
    logger.info("\n" + "="*60)
    logger.info("Curation complete!")
    logger.info("="*60)
    logger.info("""
Next steps:
1. Review curated genes in benchmarks/evidence_curation/
2. Map genes to GWAS loci using map_genes_to_loci.py
3. Rebuild Task A gold standard with expanded evidence
""")


if __name__ == '__main__':
    main()
