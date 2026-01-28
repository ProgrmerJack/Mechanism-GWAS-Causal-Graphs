"""
Download Official Baseline Scores and Benchmarking Data
==========================================================

Downloads:
1. Open Targets Genetics L2G scores (official ML predictions)
2. Open Targets gold standard curation
3. cS2G scores and code (Gazal et al. 2022)
4. FLAMES code/data (Schipper et al. 2025)
5. Fine-mapped credible sets for benchmarking

Critical: Use OFFICIAL published scores, not reimplementations.
"""

import requests
import pandas as pd
import json
import gzip
from pathlib import Path
from typing import Dict, List, Optional
import logging
import time

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class OfficialBaselineDownloader:
    """Download official baseline scores and benchmark data."""
    
    def __init__(self, data_dir: str = "data/external"):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        self.l2g_dir = self.data_dir / "open_targets"
        self.cs2g_dir = self.data_dir / "cs2g"
        self.flames_dir = self.data_dir / "flames"
        self.benchmark_dir = self.data_dir / "benchmark_gold_standards"
        
        for dir_path in [self.l2g_dir, self.cs2g_dir, self.flames_dir, self.benchmark_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
    
    def download_open_targets_l2g_scores(
        self,
        trait_efo_ids: Optional[List[str]] = None
    ) -> None:
        """
        Download official Open Targets L2G scores via GraphQL API.
        
        Uses official predictions from Mountjoy et al. (2021) Nat Genet.
        
        Args:
            trait_efo_ids: List of EFO trait IDs. If None, downloads lipid traits.
        """
        logger.info("Downloading Open Targets L2G scores...")
        
        # Default to lipid/cardiovascular traits
        if trait_efo_ids is None:
            trait_efo_ids = [
                "EFO_0004611",  # LDL cholesterol
                "EFO_0004612",  # HDL cholesterol
                "EFO_0004465",  # Total cholesterol
                "EFO_0004574",  # Triglycerides
                "EFO_0000378",  # Coronary artery disease
            ]
        
        # GraphQL endpoint
        url = "https://api.genetics.opentargets.org/graphql"
        
        all_l2g_scores = []
        
        for efo_id in trait_efo_ids:
            logger.info(f"  Fetching L2G scores for {efo_id}...")
            
            # GraphQL query for L2G scores
            query = """
            query l2gQuery($studyId: String!) {
              studyLocus2GeneTable(studyId: $studyId) {
                rows {
                  gene {
                    id
                    symbol
                  }
                  studyLocus {
                    variantId
                    chromosome
                    position
                  }
                  yProbaModel
                  yProbaDistance
                  yProbaInteraction
                  yProbaMolecularQTL
                  yProbaPathogenicity
                  hasColoc
                  distanceToLocus
                }
              }
            }
            """
            
            variables = {"studyId": efo_id}
            
            try:
                response = requests.post(
                    url,
                    json={"query": query, "variables": variables},
                    timeout=30
                )
                response.raise_for_status()
                
                data = response.json()
                
                if "data" in data and data["data"]["studyLocus2GeneTable"]:
                    rows = data["data"]["studyLocus2GeneTable"]["rows"]
                    
                    for row in rows:
                        all_l2g_scores.append({
                            "trait_efo": efo_id,
                            "gene_id": row["gene"]["id"],
                            "gene_symbol": row["gene"]["symbol"],
                            "variant_id": row["studyLocus"]["variantId"],
                            "chromosome": row["studyLocus"]["chromosome"],
                            "position": row["studyLocus"]["position"],
                            "l2g_score": row["yProbaModel"],
                            "distance_score": row["yProbaDistance"],
                            "interaction_score": row["yProbaInteraction"],
                            "qtl_score": row["yProbaMolecularQTL"],
                            "pathogenicity_score": row["yProbaPathogenicity"],
                            "has_coloc": row["hasColoc"],
                            "distance_to_locus": row["distanceToLocus"],
                        })
                    
                    logger.info(f"    Retrieved {len(rows)} L2G predictions")
                
                time.sleep(1)  # Rate limiting
                
            except Exception as e:
                logger.error(f"    Failed to download {efo_id}: {e}")
        
        if all_l2g_scores:
            df = pd.DataFrame(all_l2g_scores)
            output_file = self.l2g_dir / "open_targets_l2g_scores.tsv"
            df.to_csv(output_file, sep="\t", index=False)
            logger.info(f"✓ Saved {len(df)} L2G scores to {output_file}")
            
            # Save summary
            summary_file = self.l2g_dir / "l2g_scores_summary.txt"
            with open(summary_file, "w") as f:
                f.write("Open Targets L2G Scores Summary\n")
                f.write("=" * 60 + "\n\n")
                f.write(f"Source: Open Targets Genetics GraphQL API\n")
                f.write(f"Reference: Mountjoy et al. (2021) Nat Genet 53:1527-1533\n")
                f.write(f"DOI: 10.1038/s41588-021-00945-5\n\n")
                f.write(f"Total predictions: {len(df)}\n")
                f.write(f"Unique genes: {df['gene_symbol'].nunique()}\n")
                f.write(f"Unique loci: {df['variant_id'].nunique()}\n")
                f.write(f"Traits: {', '.join(trait_efo_ids)}\n\n")
                f.write("These are OFFICIAL L2G scores from the published model.\n")
                f.write("Do NOT claim to have 'implemented L2G' - use these scores directly.\n")
    
    def download_open_targets_gold_standard(self) -> None:
        """Download Open Targets curated gene-disease associations."""
        logger.info("Downloading Open Targets gold standard curation...")
        
        url = "https://api.genetics.opentargets.org/graphql"
        
        # Query for high-confidence gene-disease associations
        query = """
        query goldStandard {
          diseases(page: {index: 0, size: 1000}) {
            rows {
              id
              name
              knownDrugs {
                rows {
                  drugId
                  drugName
                  targetId
                  targetName
                  phase
                }
              }
            }
          }
        }
        """
        
        try:
            response = requests.post(url, json={"query": query}, timeout=60)
            response.raise_for_status()
            
            data = response.json()
            
            gold_standard_genes = []
            
            if "data" in data and "diseases" in data["data"]:
                diseases = data["data"]["diseases"]["rows"]
                
                for disease in diseases:
                    if disease["knownDrugs"] and disease["knownDrugs"]["rows"]:
                        for drug in disease["knownDrugs"]["rows"]:
                            if drug["phase"] >= 3:  # Phase 3+ approved drugs
                                gold_standard_genes.append({
                                    "disease_id": disease["id"],
                                    "disease_name": disease["name"],
                                    "gene_id": drug["targetId"],
                                    "gene_symbol": drug["targetName"],
                                    "drug_id": drug["drugId"],
                                    "drug_name": drug["drugName"],
                                    "clinical_phase": drug["phase"],
                                    "evidence_type": "approved_drug_target",
                                })
            
            if gold_standard_genes:
                df = pd.DataFrame(gold_standard_genes)
                output_file = self.benchmark_dir / "open_targets_gold_standard.tsv"
                df.to_csv(output_file, sep="\t", index=False)
                logger.info(f"✓ Saved {len(df)} gold standard genes to {output_file}")
        
        except Exception as e:
            logger.error(f"Failed to download gold standard: {e}")
    
    def download_cs2g_resources(self) -> None:
        """
        Download cS2G code and documentation.
        
        Note: cS2G requires heritability enrichment calculations that need
        custom GWAS data. We document the official source but don't reimplement.
        """
        logger.info("Documenting cS2G resources...")
        
        info_file = self.cs2g_dir / "CS2G_INFO.md"
        with open(info_file, "w") as f:
            f.write("# cS2G (Cross-ancestry SNP-to-gene) Resources\n\n")
            f.write("**Reference:** Gazal et al. (2022) Nature Genetics 54:707-717\n")
            f.write("**DOI:** 10.1038/s41588-022-01087-y\n\n")
            f.write("## Official Code\n")
            f.write("- Repository: https://alkesgroup.broadinstitute.org/cS2G/code\n")
            f.write("- GitHub: https://github.com/omerwe/cS2G\n\n")
            f.write("## Method Summary\n")
            f.write("cS2G combines 7 SNP-to-gene linking strategies weighted by\n")
            f.write("heritability enrichment:\n")
            f.write("1. ABC enhancer-gene links\n")
            f.write("2. eQTL colocalization\n")
            f.write("3. Distance to TSS\n")
            f.write("4. Coding variants\n")
            f.write("5. Chromatin interactions (PCHi-C, HiC)\n")
            f.write("6. ATAC-seq peaks\n")
            f.write("7. Evolutionary constraint\n\n")
            f.write("## Implementation Notes\n")
            f.write("cS2G requires:\n")
            f.write("- GWAS summary statistics\n")
            f.write("- Fine-mapped credible sets\n")
            f.write("- Heritability enrichment calculations (LD score regression)\n")
            f.write("- Pre-computed SNP-to-gene annotations\n\n")
            f.write("**DO NOT claim 'we implemented cS2G' unless:**\n")
            f.write("1. You ran their official code, OR\n")
            f.write("2. You exactly reproduced their heritability weighting scheme\n\n")
            f.write("**Alternative phrasing:**\n")
            f.write("- 'cS2G-inspired baseline using ABC + eQTL + distance'\n")
            f.write("- 'ABC+eQTL composite score (similar to cS2G strategy)'\n")
            f.write("- Benchmark against published cS2G results (cite paper)\n")
        
        logger.info(f"✓ Documented cS2G at {info_file}")
    
    def download_flames_resources(self) -> None:
        """
        Document FLAMES (Schipper et al. 2025) resources.
        
        FLAMES is a complex pipeline requiring PoPS integration and XGBoost training.
        """
        logger.info("Documenting FLAMES resources...")
        
        info_file = self.flames_dir / "FLAMES_INFO.md"
        with open(info_file, "w") as f:
            f.write("# FLAMES (Fine-mapped Locus Assessment Model of Effector genes)\n\n")
            f.write("**Reference:** Schipper et al. (2025) Nature Genetics 57(11):2930\n")
            f.write("**DOI:** 10.1038/s41588-025-02416-7\n")
            f.write("**Preprint:** bioRxiv 2023.12.23.23300360\n\n")
            f.write("## Official Code\n")
            f.write("- GitHub: https://github.com/Marijn-Schipper/FLAMES\n")
            f.write("- Analysis code: https://github.com/Marijn-Schipper/FLAMES_paper_analyses\n")
            f.write("- Reference data: https://zenodo.org/uploads/10409723\n\n")
            f.write("## Method Summary\n")
            f.write("FLAMES combines:\n")
            f.write("1. XGBoost classifier trained on curated gene-locus pairs\n")
            f.write("   - Features: ABC, eQTL, distance, coding, PIPs\n")
            f.write("2. PoPS convergence scores (gene-level pathway enrichment)\n")
            f.write("3. Integration: 0.725 * XGBoost + 0.275 * PoPS\n\n")
            f.write("## Implementation Requirements\n")
            f.write("FLAMES requires:\n")
            f.write("- Fine-mapped credible sets with PIPs\n")
            f.write("- ABC enhancer-gene predictions\n")
            f.write("- eQTL colocalization scores\n")
            f.write("- PoPS gene scores (pathway-naïve version)\n")
            f.write("- XGBoost >=1.6.0\n")
            f.write("- Curated training data (gene-locus pairs)\n\n")
            f.write("**DO NOT claim 'we implemented FLAMES' unless:**\n")
            f.write("1. You ran their official code from GitHub, OR\n")
            f.write("2. You trained XGBoost with their exact feature set AND integrated with PoPS\n\n")
            f.write("**Alternative phrasing:**\n")
            f.write("- 'FLAMES-inspired XGBoost baseline'\n")
            f.write("- 'ABC+eQTL+PIP gradient boosting (similar to FLAMES)'\n")
            f.write("- Benchmark against published FLAMES results (cite paper)\n")
        
        logger.info(f"✓ Documented FLAMES at {info_file}")
    
    def download_finemapping_gold_standard(self) -> None:
        """
        Download fine-mapped coding variants as unambiguous gold standard.
        
        These are loci where the credible set contains a missense/loss-of-function
        variant, making the causal gene unambiguous.
        """
        logger.info("Creating fine-mapped coding variant gold standard...")
        
        # This would typically come from a curated database
        # For now, create template showing structure
        template_file = self.benchmark_dir / "finemapped_coding_gold_standard_TEMPLATE.tsv"
        
        template_df = pd.DataFrame({
            "locus_id": ["LDLR_19:11200038", "PCSK9_1:55039974", "APOB_2:21001429"],
            "gene_symbol": ["LDLR", "PCSK9", "APOB"],
            "gene_id": ["ENSG00000130203", "ENSG00000169174", "ENSG00000084674"],
            "lead_variant": ["rs688", "rs505151", "rs1367117"],
            "chromosome": ["19", "1", "2"],
            "position": [11200038, 55039974, 21001429],
            "trait": ["LDL_cholesterol", "LDL_cholesterol", "LDL_cholesterol"],
            "coding_variant_in_credible_set": ["rs688_missense", "rs505151_missense", "rs1367117_missense"],
            "pip": [0.99, 0.97, 0.95],
            "evidence_type": ["finemapped_coding", "finemapped_coding", "finemapped_coding"],
            "confidence": [1.0, 1.0, 1.0],
        })
        
        template_df.to_csv(template_file, sep="\t", index=False)
        logger.info(f"✓ Created gold standard template at {template_file}")
        logger.info("  NOTE: Fill this with real fine-mapped credible sets")
    
    def run_all_downloads(self) -> None:
        """Execute all downloads."""
        logger.info("=" * 70)
        logger.info("DOWNLOADING OFFICIAL BASELINE SCORES AND GOLD STANDARDS")
        logger.info("=" * 70)
        
        # Download Open Targets L2G scores (official)
        self.download_open_targets_l2g_scores()
        
        # Download Open Targets gold standard
        self.download_open_targets_gold_standard()
        
        # Document cS2G and FLAMES (complex pipelines)
        self.download_cs2g_resources()
        self.download_flames_resources()
        
        # Create gold standard templates
        self.download_finemapping_gold_standard()
        
        logger.info("\n" + "=" * 70)
        logger.info("DOWNLOAD COMPLETE")
        logger.info("=" * 70)
        logger.info("\nNext steps:")
        logger.info("1. Use official L2G scores (data/external/open_targets/)")
        logger.info("2. Rename any 'approximations' in baseline code")
        logger.info("3. Add baseline provenance table to manuscript")
        logger.info("4. Run benchmarking with proper gold standards")


if __name__ == "__main__":
    downloader = OfficialBaselineDownloader()
    downloader.run_all_downloads()
