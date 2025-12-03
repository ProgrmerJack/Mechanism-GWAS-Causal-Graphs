"""
Tests for Mechanism Graph Module
"""

import pytest
import numpy as np

from src.mechanism_graph.nodes import (
    VariantNode,
    CCRENode,
    GeneNode,
    TissueNode,
    TraitNode,
)
from src.mechanism_graph.edges import (
    EdgeProbability,
    VariantToCCREEdge,
    CCREToGeneEdge,
    GeneToTissueEdge,
)
from src.mechanism_graph.graph import MechanismGraph


class TestNodes:
    """Tests for node types."""
    
    def test_variant_node(self):
        """Test VariantNode creation."""
        node = VariantNode(
            rsid="rs12740374",
            chrom="1",
            pos=109817590,
            ref="G",
            alt="T",
            pip=0.94,
            cs_id=1,
        )
        
        assert node.rsid == "rs12740374"
        assert node.pip == 0.94
        assert node.node_type == "variant"
    
    def test_ccre_node(self):
        """Test CCRENode creation."""
        node = CCRENode(
            ccre_id="EH38E2387531",
            chrom="1",
            start=109817000,
            end=109818000,
            element_type="enhancer",
            tissue_activities={"liver": 0.9, "heart": 0.3},
        )
        
        assert node.element_type == "enhancer"
        assert node.tissue_activities["liver"] == 0.9
    
    def test_gene_node(self):
        """Test GeneNode creation."""
        node = GeneNode(
            gene_id="ENSG00000134222",
            symbol="SORT1",
            chrom="1",
            start=109702000,
            end=109815000,
            strand="+",
            biotype="protein_coding",
            expression={"liver": 50.5, "heart": 10.2},
        )
        
        assert node.symbol == "SORT1"
        assert node.expression["liver"] == 50.5
    
    def test_tissue_node(self):
        """Test TissueNode creation."""
        node = TissueNode(
            tissue_id="GTEX_LIVER",
            name="Liver",
            ontology_id="UBERON:0002107",
            relevance_prior=0.9,
        )
        
        assert node.name == "Liver"
        assert node.relevance_prior == 0.9
    
    def test_trait_node(self):
        """Test TraitNode creation."""
        node = TraitNode(
            trait_id="ldl_cholesterol",
            name="LDL Cholesterol",
            category="lipid",
        )
        
        assert node.name == "LDL Cholesterol"


class TestEdges:
    """Tests for edge probability computation."""
    
    def test_variant_to_ccre_overlap(self):
        """Test variant-cCRE overlap edge."""
        edge = VariantToCCREEdge()
        
        # Perfect overlap
        prob = edge.compute(
            variant_pos=109817500,
            ccre_start=109817000,
            ccre_end=109818000,
            variant_pip=0.9,
        )
        
        assert prob.probability > 0.8
        assert prob.source == "overlap"
    
    def test_variant_to_ccre_distance(self):
        """Test variant-cCRE distance decay."""
        edge = VariantToCCREEdge()
        
        # Far from cCRE
        prob = edge.compute(
            variant_pos=109900000,
            ccre_start=109817000,
            ccre_end=109818000,
            variant_pip=0.9,
        )
        
        assert prob.probability < 0.5
    
    def test_ccre_to_gene_distance(self):
        """Test cCRE-gene distance edge."""
        edge = CCREToGeneEdge()
        
        # Close to TSS
        prob = edge.compute(
            ccre_center=109700000,
            gene_tss=109702000,
            abc_score=None,
        )
        
        assert prob.probability > 0.5
    
    def test_ccre_to_gene_abc(self):
        """Test cCRE-gene ABC score edge."""
        edge = CCREToGeneEdge()
        
        # With ABC score
        prob = edge.compute(
            ccre_center=109700000,
            gene_tss=109702000,
            abc_score=0.8,
        )
        
        assert prob.probability > 0.7
        assert "abc" in prob.source.lower()
    
    def test_gene_to_tissue(self):
        """Test gene-tissue colocalization edge."""
        edge = GeneToTissueEdge()
        
        prob = edge.compute(
            coloc_h4=0.92,
            n_variants=5,
        )
        
        assert prob.probability == 0.92


class TestMechanismGraph:
    """Tests for graph construction and inference."""
    
    def test_graph_creation(self):
        """Test basic graph creation."""
        graph = MechanismGraph(
            locus_id="chr1_109817590",
            trait="ldl_cholesterol",
        )
        
        assert graph.locus_id == "chr1_109817590"
        assert len(graph.nodes) == 0
    
    def test_add_node(self):
        """Test adding nodes to graph."""
        graph = MechanismGraph(
            locus_id="chr1_109817590",
            trait="ldl_cholesterol",
        )
        
        variant = VariantNode(
            rsid="rs12740374",
            chrom="1",
            pos=109817590,
            ref="G",
            alt="T",
            pip=0.94,
        )
        
        graph.add_node(variant)
        
        assert len(graph.nodes) == 1
        assert graph.nodes["rs12740374"] == variant
    
    def test_add_edge(self):
        """Test adding edges to graph."""
        graph = MechanismGraph(
            locus_id="chr1_109817590",
            trait="ldl_cholesterol",
        )
        
        variant = VariantNode(
            rsid="rs12740374",
            chrom="1",
            pos=109817590,
            ref="G",
            alt="T",
            pip=0.94,
        )
        
        ccre = CCRENode(
            ccre_id="EH38E2387531",
            chrom="1",
            start=109817000,
            end=109818000,
            element_type="enhancer",
        )
        
        graph.add_node(variant)
        graph.add_node(ccre)
        
        edge_prob = EdgeProbability(
            source_id="rs12740374",
            target_id="EH38E2387531",
            probability=0.9,
            confidence=0.85,
            source="overlap",
        )
        
        graph.add_edge(edge_prob)
        
        assert len(graph.edges) == 1
    
    def test_gene_score_computation(self):
        """Test gene score computation from graph."""
        graph = MechanismGraph(
            locus_id="chr1_109817590",
            trait="ldl_cholesterol",
        )
        
        # Add nodes
        variant = VariantNode(
            rsid="rs12740374",
            chrom="1",
            pos=109817590,
            ref="G",
            alt="T",
            pip=0.94,
        )
        
        ccre = CCRENode(
            ccre_id="EH38E2387531",
            chrom="1",
            start=109817000,
            end=109818000,
            element_type="enhancer",
        )
        
        gene = GeneNode(
            gene_id="ENSG00000134222",
            symbol="SORT1",
            chrom="1",
            start=109702000,
            end=109815000,
            strand="+",
            biotype="protein_coding",
        )
        
        tissue = TissueNode(
            tissue_id="liver",
            name="Liver",
            relevance_prior=0.9,
        )
        
        trait = TraitNode(
            trait_id="ldl_cholesterol",
            name="LDL Cholesterol",
        )
        
        for node in [variant, ccre, gene, tissue, trait]:
            graph.add_node(node)
        
        # Add edges
        edges = [
            EdgeProbability("rs12740374", "EH38E2387531", 0.9, 0.8, "overlap"),
            EdgeProbability("EH38E2387531", "ENSG00000134222", 0.8, 0.7, "distance"),
            EdgeProbability("ENSG00000134222", "liver", 0.92, 0.85, "coloc"),
            EdgeProbability("liver", "ldl_cholesterol", 0.9, 0.9, "prior"),
        ]
        
        for edge in edges:
            graph.add_edge(edge)
        
        # Compute gene score
        scores = graph.compute_gene_scores()
        
        assert "SORT1" in scores or "ENSG00000134222" in scores
    
    def test_graph_serialization(self):
        """Test graph JSON serialization."""
        graph = MechanismGraph(
            locus_id="chr1_109817590",
            trait="ldl_cholesterol",
        )
        
        variant = VariantNode(
            rsid="rs12740374",
            chrom="1",
            pos=109817590,
            ref="G",
            alt="T",
            pip=0.94,
        )
        
        graph.add_node(variant)
        
        json_data = graph.to_json()
        
        assert "locus_id" in json_data
        assert "nodes" in json_data


class TestCalibration:
    """Tests for calibration module."""
    
    def test_reliability_curve(self):
        """Test reliability curve computation."""
        from src.calibration.calibration import CalibrationAnalysis
        
        np.random.seed(42)
        
        # Generate well-calibrated predictions
        n = 1000
        y_prob = np.random.beta(2, 5, n)
        y_true = (np.random.random(n) < y_prob).astype(int)
        
        analyzer = CalibrationAnalysis(n_bins=10)
        curve = analyzer.compute_reliability_curve(y_true, y_prob)
        
        assert "mean_predicted_prob" in curve
        assert "fraction_positives" in curve
        assert len(curve["mean_predicted_prob"]) <= 10
    
    def test_ece(self):
        """Test Expected Calibration Error."""
        from src.calibration.calibration import CalibrationAnalysis
        
        np.random.seed(42)
        
        # Perfect calibration
        n = 1000
        y_prob = np.random.random(n)
        y_true = (np.random.random(n) < y_prob).astype(int)
        
        analyzer = CalibrationAnalysis(n_bins=10)
        ece = analyzer.expected_calibration_error(y_true, y_prob)
        
        # Should be low for well-calibrated predictions
        assert ece < 0.2
    
    def test_brier_decomposition(self):
        """Test Brier score decomposition."""
        from src.calibration.calibration import CalibrationAnalysis
        
        np.random.seed(42)
        
        n = 1000
        y_prob = np.random.random(n)
        y_true = (np.random.random(n) < y_prob).astype(int)
        
        analyzer = CalibrationAnalysis()
        decomp = analyzer.brier_score_decomposition(y_true, y_prob)
        
        assert "brier_score" in decomp
        assert "reliability" in decomp
        assert "resolution" in decomp
        assert "uncertainty" in decomp


class TestBenchmarks:
    """Tests for benchmark module."""
    
    def test_gold_standard_loading(self):
        """Test gold standard gene loading."""
        from src.calibration.benchmarks import GoldStandardGenes
        
        gs = GoldStandardGenes("ldl_cholesterol")
        
        assert "LDLR" in gs.genes
        assert "PCSK9" in gs.genes
        assert gs.is_gold_standard("LDLR")
    
    def test_benchmark_evaluation(self):
        """Test benchmark evaluation."""
        from src.calibration.benchmarks import BenchmarkSuite
        
        suite = BenchmarkSuite("ldl_cholesterol")
        
        # Create mock rankings
        rankings = [
            {"gene": "LDLR", "score": 0.95},
            {"gene": "PCSK9", "score": 0.90},
            {"gene": "APOB", "score": 0.85},
            {"gene": "FAKE1", "score": 0.80},
            {"gene": "FAKE2", "score": 0.75},
        ] + [{"gene": f"GENE{i}", "score": 0.5 - i*0.01} for i in range(95)]
        
        results = suite.evaluate_ranking(rankings)
        
        assert "recall@10" in results
        assert results["recall@10"] > 0  # Should find some gold standard genes


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
