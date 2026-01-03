#!/usr/bin/env python3
"""
Comprehensive Data Validation for Mechanism-GWAS-Causal-Graphs Project

This script performs full-scale statistical analysis on ALL available data
to validate manuscript claims with deep statistical rigor.
"""

import json
import sys
import os
from pathlib import Path
from collections import defaultdict
import statistics
import math

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))


def load_tsv(filepath):
    """Load TSV file into list of dictionaries."""
    records = []
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        if not lines:
            return records
        headers = lines[0].strip().split('\t')
        for line in lines[1:]:
            if line.strip():
                values = line.strip().split('\t')
                record = {}
                for i, h in enumerate(headers):
                    record[h] = values[i] if i < len(values) else ''
                records.append(record)
    return records


def load_json(filepath):
    """Load JSON file."""
    with open(filepath, 'r', encoding='utf-8') as f:
        return json.load(f)


def load_yaml(filepath):
    """Load YAML file as simple dict (basic parser)."""
    data = {}
    current_key = None
    current_section = None
    
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#') or not line.strip():
                continue
            
            # Simple key-value parsing
            if ':' in line and not line.startswith(' '):
                parts = line.split(':', 1)
                key = parts[0].strip()
                value = parts[1].strip() if len(parts) > 1 else ''
                if value:
                    # Try to parse numbers
                    try:
                        if '.' in value:
                            data[key] = float(value)
                        else:
                            data[key] = int(value)
                    except ValueError:
                        data[key] = value.strip('"').strip("'")
                else:
                    current_section = key
                    data[current_section] = {}
            elif current_section and ':' in line:
                parts = line.strip().split(':', 1)
                key = parts[0].strip()
                value = parts[1].strip() if len(parts) > 1 else ''
                try:
                    if '.' in value:
                        data[current_section][key] = float(value)
                    elif value.isdigit():
                        data[current_section][key] = int(value)
                    else:
                        data[current_section][key] = value.strip('"').strip("'")
                except:
                    pass
    return data


def compute_calibration_statistics(reliability_data):
    """Compute comprehensive calibration statistics from reliability diagram data."""
    results = {}
    
    # Group by module
    modules = defaultdict(list)
    for record in reliability_data:
        module = record.get('module', 'unknown')
        modules[module].append(record)
    
    for module, bins in modules.items():
        # Calculate ECE (Expected Calibration Error)
        total_samples = sum(int(b.get('n_samples', 0)) for b in bins)
        ece = 0
        mce = 0  # Maximum Calibration Error
        brier_decomp = {'reliability': 0, 'resolution': 0, 'uncertainty': 0}
        
        predictions = []
        outcomes = []
        
        for b in bins:
            n = int(b.get('n_samples', 0))
            pred = float(b.get('mean_predicted', 0))
            obs = float(b.get('observed_frequency', 0))
            
            if total_samples > 0:
                weight = n / total_samples
                calibration_error = abs(pred - obs)
                ece += weight * calibration_error
                mce = max(mce, calibration_error)
                
                # For Brier decomposition
                predictions.extend([pred] * n)
                outcomes.extend([obs] * n)
        
        results[module] = {
            'ECE': round(ece, 4),
            'MCE': round(mce, 4),
            'n_bins': len(bins),
            'total_samples': total_samples,
            'calibration_status': 'EXCELLENT' if ece < 0.05 else 'GOOD' if ece < 0.1 else 'MODERATE'
        }
    
    return results


def analyze_benchmark_genes(tier1_data, tier2_data):
    """Analyze gold standard genes and drug targets comprehensively."""
    results = {
        'tier1_analysis': {},
        'tier2_analysis': {},
        'combined_analysis': {}
    }
    
    # Tier 1 analysis
    traits = defaultdict(list)
    evidence_types = defaultdict(int)
    validation_counts = {'anti_leak': 0, 'training_excluded': 0}
    
    for gene in tier1_data:
        trait = gene.get('trait', 'unknown')
        traits[trait].append(gene)
        
        evidence = gene.get('evidence_type', '')
        for ev_type in evidence.split(';'):
            evidence_types[ev_type.strip()] += 1
        
        if gene.get('anti_leak_verified', '').upper() == 'TRUE':
            validation_counts['anti_leak'] += 1
        if gene.get('training_set_exclusion_verified', '').upper() == 'TRUE':
            validation_counts['training_excluded'] += 1
    
    results['tier1_analysis'] = {
        'total_genes': len(tier1_data),
        'traits_covered': list(traits.keys()),
        'genes_per_trait': {t: len(g) for t, g in traits.items()},
        'evidence_type_distribution': dict(evidence_types),
        'anti_leak_verified_pct': round(validation_counts['anti_leak'] / len(tier1_data) * 100, 1) if tier1_data else 0,
        'training_excluded_pct': round(validation_counts['training_excluded'] / len(tier1_data) * 100, 1) if tier1_data else 0
    }
    
    # Tier 2 analysis (drug targets)
    drugs_per_gene = defaultdict(list)
    approval_status = defaultdict(int)
    mechanisms = defaultdict(int)
    
    for target in tier2_data:
        gene = target.get('gene_symbol', 'unknown')
        drug = target.get('drug_name', 'unknown')
        drugs_per_gene[gene].append(drug)
        
        status = target.get('approval_status', 'Unknown')
        approval_status[status] += 1
        
        mechanism = target.get('mechanism_of_action', 'unknown')
        mechanisms[mechanism] += 1
    
    results['tier2_analysis'] = {
        'total_drugs': len(tier2_data),
        'unique_genes': len(drugs_per_gene),
        'drugs_per_gene': {g: len(d) for g, d in drugs_per_gene.items()},
        'approval_status_distribution': dict(approval_status),
        'mechanism_distribution': dict(mechanisms)
    }
    
    # Combined statistics
    tier1_genes = set(g['gene_symbol'] for g in tier1_data if 'gene_symbol' in g)
    tier2_genes = set(t['gene_symbol'] for t in tier2_data if 'gene_symbol' in t)
    
    results['combined_analysis'] = {
        'tier1_genes': len(tier1_genes),
        'tier2_genes': len(tier2_genes),
        'overlap': len(tier1_genes & tier2_genes),
        'tier1_only': len(tier1_genes - tier2_genes),
        'tier2_only': len(tier2_genes - tier1_genes),
        'union': len(tier1_genes | tier2_genes)
    }
    
    return results


def analyze_replication_data(replication_data, summary_yaml):
    """Comprehensive analysis of eQTL replication data."""
    results = {
        'gene_level': {},
        'tissue_level': {},
        'effect_size_analysis': {},
        'manuscript_claim_validation': {}
    }
    
    # Gene-level analysis
    genes_replicated = 0
    genes_total = len(replication_data)
    pp_h4_gtex = []
    pp_h4_repl = []
    effect_correlations = []
    direction_concordant = 0
    
    tissue_stats = defaultdict(lambda: {'replicated': 0, 'total': 0})
    
    for gene in replication_data:
        replicated = gene.get('replicated', '').upper() == 'TRUE'
        if replicated:
            genes_replicated += 1
        
        tissue = gene.get('gtex_tissue', 'unknown')
        tissue_stats[tissue]['total'] += 1
        if replicated:
            tissue_stats[tissue]['replicated'] += 1
        
        try:
            gtex_h4 = float(gene.get('gtex_pp_h4', 0))
            repl_h4 = float(gene.get('replication_pp_h4', 0))
            pp_h4_gtex.append(gtex_h4)
            pp_h4_repl.append(repl_h4)
        except:
            pass
        
        try:
            corr = float(gene.get('effect_correlation', 0))
            effect_correlations.append(corr)
        except:
            pass
        
        if gene.get('direction_concordant', '').upper() == 'TRUE':
            direction_concordant += 1
    
    # Compute statistics
    replication_rate = genes_replicated / genes_total if genes_total > 0 else 0
    
    results['gene_level'] = {
        'total_genes': genes_total,
        'genes_replicated': genes_replicated,
        'replication_rate': round(replication_rate * 100, 1),
        'direction_concordance_rate': round(direction_concordant / genes_total * 100, 1) if genes_total > 0 else 0
    }
    
    results['tissue_level'] = {
        tissue: {
            'replicated': stats['replicated'],
            'total': stats['total'],
            'rate': round(stats['replicated'] / stats['total'] * 100, 1) if stats['total'] > 0 else 0
        }
        for tissue, stats in tissue_stats.items()
    }
    
    if effect_correlations:
        results['effect_size_analysis'] = {
            'mean_correlation': round(statistics.mean(effect_correlations), 3),
            'median_correlation': round(statistics.median(effect_correlations), 3),
            'min_correlation': round(min(effect_correlations), 3),
            'max_correlation': round(max(effect_correlations), 3),
            'std_correlation': round(statistics.stdev(effect_correlations), 3) if len(effect_correlations) > 1 else 0
        }
    
    if pp_h4_gtex and pp_h4_repl:
        # Compute correlation between GTEx and replication PP.H4
        n = len(pp_h4_gtex)
        mean_gtex = sum(pp_h4_gtex) / n
        mean_repl = sum(pp_h4_repl) / n
        
        numerator = sum((g - mean_gtex) * (r - mean_repl) for g, r in zip(pp_h4_gtex, pp_h4_repl))
        denom_gtex = math.sqrt(sum((g - mean_gtex)**2 for g in pp_h4_gtex))
        denom_repl = math.sqrt(sum((r - mean_repl)**2 for r in pp_h4_repl))
        
        correlation = numerator / (denom_gtex * denom_repl) if denom_gtex * denom_repl > 0 else 0
        
        results['effect_size_analysis']['pp_h4_correlation'] = round(correlation, 3)
    
    # Validate manuscript claims
    results['manuscript_claim_validation'] = {
        'claimed_replication_rate_78pct': True,
        'actual_replication_rate': round(replication_rate * 100, 1),
        'exceeds_claim': replication_rate >= 0.78,
        'status': 'VALIDATED' if replication_rate >= 0.78 else 'BELOW_CLAIM'
    }
    
    return results


def analyze_mechanism_graphs(mechanism_graphs):
    """Analyze mechanism graph structures comprehensively."""
    results = {
        'overall': {},
        'per_locus': {}
    }
    
    total_paths = 0
    total_genes = 0
    path_probabilities = []
    coloc_pp_h4 = []
    replicated_count = 0
    
    for locus_name, graph in mechanism_graphs.items():
        locus_results = {
            'locus_id': graph.get('locus_id', 'unknown'),
            'trait': graph.get('trait', 'unknown'),
            'n_credible_variants': len(graph.get('credible_variants', [])),
            'n_mechanism_paths': len(graph.get('mechanism_paths', [])),
            'n_genes': len(graph.get('gene_probabilities', {}))
        }
        
        # Analyze credible variants
        variants = graph.get('credible_variants', [])
        max_pip = max((float(v.get('pip', 0)) for v in variants), default=0)
        locus_results['max_variant_pip'] = round(max_pip, 3)
        
        # Analyze mechanism paths
        paths = graph.get('mechanism_paths', [])
        for path in paths:
            total_paths += 1
            prob = float(path.get('path_probability', 0))
            path_probabilities.append(prob)
            
            coloc = path.get('colocalization', {})
            pp_h4 = float(coloc.get('pp_h4', 0))
            if pp_h4 > 0:
                coloc_pp_h4.append(pp_h4)
            
            replication = path.get('replication', {})
            if replication.get('replicated', False):
                replicated_count += 1
        
        # Analyze gene probabilities
        genes = graph.get('gene_probabilities', {})
        total_genes += len(genes)
        if genes:
            top_gene = max(genes.items(), key=lambda x: x[1].get('probability', 0))
            locus_results['top_gene'] = top_gene[0]
            locus_results['top_gene_probability'] = round(top_gene[1].get('probability', 0), 3)
        
        # Validation status
        validation = graph.get('validation', {})
        locus_results['is_gold_standard'] = validation.get('tier1_gold_standard', False)
        locus_results['prediction_correct'] = validation.get('prediction_correct', False)
        
        results['per_locus'][locus_name] = locus_results
    
    # Overall statistics
    results['overall'] = {
        'n_loci': len(mechanism_graphs),
        'total_mechanism_paths': total_paths,
        'total_genes_analyzed': total_genes,
        'avg_path_probability': round(statistics.mean(path_probabilities), 3) if path_probabilities else 0,
        'avg_coloc_pp_h4': round(statistics.mean(coloc_pp_h4), 3) if coloc_pp_h4 else 0,
        'paths_replicated': replicated_count,
        'replication_rate': round(replicated_count / total_paths * 100, 1) if total_paths > 0 else 0
    }
    
    return results


def analyze_locus_summary(locus_data):
    """Comprehensive analysis of locus summary data."""
    results = {
        'trait_distribution': {},
        'evidence_tier_distribution': {},
        'validation_status_distribution': {},
        'mechanism_type_distribution': {},
        'statistical_summary': {}
    }
    
    traits = defaultdict(int)
    evidence_tiers = defaultdict(int)
    validation_status = defaultdict(int)
    mechanism_types = defaultdict(int)
    
    path_probs = []
    coloc_pp_h4s = []
    max_pips = []
    n_genes_tested = []
    
    for locus in locus_data:
        trait = locus.get('trait', 'unknown')
        traits[trait] += 1
        
        tier = locus.get('evidence_tier', 'unknown')
        evidence_tiers[tier] += 1
        
        status = locus.get('validation_status', 'unknown')
        validation_status[status] += 1
        
        mechanism = locus.get('mechanism_type', 'unknown')
        mechanism_types[mechanism] += 1
        
        try:
            path_probs.append(float(locus.get('path_probability', 0)))
        except:
            pass
        
        try:
            coloc_pp_h4s.append(float(locus.get('coloc_pp_h4', 0)))
        except:
            pass
        
        try:
            max_pips.append(float(locus.get('max_pip', 0)))
        except:
            pass
        
        try:
            n_genes_tested.append(int(locus.get('n_genes_tested', 0)))
        except:
            pass
    
    results['trait_distribution'] = dict(traits)
    results['evidence_tier_distribution'] = dict(evidence_tiers)
    results['validation_status_distribution'] = dict(validation_status)
    results['mechanism_type_distribution'] = dict(mechanism_types)
    
    if path_probs:
        results['statistical_summary']['path_probability'] = {
            'mean': round(statistics.mean(path_probs), 3),
            'median': round(statistics.median(path_probs), 3),
            'std': round(statistics.stdev(path_probs), 3) if len(path_probs) > 1 else 0,
            'min': round(min(path_probs), 3),
            'max': round(max(path_probs), 3)
        }
    
    if coloc_pp_h4s:
        results['statistical_summary']['colocalization_pp_h4'] = {
            'mean': round(statistics.mean(coloc_pp_h4s), 3),
            'median': round(statistics.median(coloc_pp_h4s), 3),
            'std': round(statistics.stdev(coloc_pp_h4s), 3) if len(coloc_pp_h4s) > 1 else 0,
            'min': round(min(coloc_pp_h4s), 3),
            'max': round(max(coloc_pp_h4s), 3)
        }
    
    if n_genes_tested:
        results['statistical_summary']['genes_per_locus'] = {
            'mean': round(statistics.mean(n_genes_tested), 1),
            'median': round(statistics.median(n_genes_tested), 1),
            'total': sum(n_genes_tested)
        }
    
    results['statistical_summary']['total_loci'] = len(locus_data)
    
    return results


def analyze_calibration_metrics(calibration_data):
    """Analyze calibration metrics comprehensively."""
    results = {
        'module_comparison': {},
        'our_method_vs_baselines': {},
        'manuscript_claims': {}
    }
    
    # Group by module
    modules = {}
    for record in calibration_data:
        module = record.get('module', 'unknown')
        metric = record.get('metric', 'unknown')
        
        if module not in modules:
            modules[module] = {}
        
        try:
            value = float(record.get('value', 0))
            ci_lower = float(record.get('ci_lower', 0))
            ci_upper = float(record.get('ci_upper', 0))
            
            modules[module][metric] = {
                'value': round(value, 3),
                'ci_95': [round(ci_lower, 3), round(ci_upper, 3)]
            }
        except:
            pass
    
    results['module_comparison'] = modules
    
    # Compare our method to baselines
    our_method = modules.get('Final_gene_probability', {})
    l2g = modules.get('Open_Targets_L2G', {})
    pops = modules.get('PoPS', {})
    magma = modules.get('MAGMA', {})
    nearest = modules.get('Nearest_gene', {})
    
    if our_method:
        comparisons = {
            'our_method': our_method,
            'vs_L2G': {},
            'vs_PoPS': {},
            'vs_MAGMA': {},
            'vs_Nearest': {}
        }
        
        for metric in ['ECE', 'Recall_at_20']:
            if metric in our_method:
                our_val = our_method[metric]['value']
                
                if metric in l2g:
                    improvement = ((l2g[metric]['value'] - our_val) / l2g[metric]['value']) * 100 if l2g[metric]['value'] != 0 else 0
                    comparisons['vs_L2G'][metric] = {
                        'ours': our_val,
                        'theirs': l2g[metric]['value'],
                        'improvement_pct': round(improvement, 1)
                    }
                
                if metric in pops:
                    improvement = ((pops[metric]['value'] - our_val) / pops[metric]['value']) * 100 if pops[metric]['value'] != 0 else 0
                    comparisons['vs_PoPS'][metric] = {
                        'ours': our_val,
                        'theirs': pops[metric]['value'],
                        'improvement_pct': round(improvement, 1)
                    }
        
        results['our_method_vs_baselines'] = comparisons
    
    # Validate manuscript claims
    claims = []
    
    # Claim 1: ECE < 0.05
    if our_method.get('ECE'):
        ece = our_method['ECE']['value']
        claims.append({
            'claim': 'ECE < 0.05 (well-calibrated)',
            'value': ece,
            'validated': ece < 0.05
        })
    
    # Claim 2: Recall@20 = 0.76
    if our_method.get('Recall_at_20'):
        recall = our_method['Recall_at_20']['value']
        claims.append({
            'claim': 'Recall@20 >= 0.76',
            'value': recall,
            'validated': recall >= 0.76
        })
    
    # Claim 3: CRISPR AUPRC = 0.71
    ccre_module = modules.get('cCRE_Gene_ABC_PCHiC', {})
    if ccre_module.get('AUPRC'):
        auprc = ccre_module['AUPRC']['value']
        claims.append({
            'claim': 'CRISPR benchmark AUPRC >= 0.71',
            'value': auprc,
            'validated': auprc >= 0.71
        })
    
    results['manuscript_claims'] = claims
    
    return results


def validate_gwas_sumstats_manifest(manifest_data, raw_data_path):
    """Validate GWAS summary statistics against manifest."""
    results = {
        'datasets_expected': [],
        'datasets_found': [],
        'validation_status': {}
    }
    
    # Parse manifest for datasets
    # This is a simplified check - just verify file existence
    gwas_dirs = ['cardiogram_2022', 'diagram_2022', 'glgc_2021', 'icbp_bp']
    
    for dir_name in gwas_dirs:
        dir_path = raw_data_path / 'gwas_sumstats' / dir_name
        if dir_path.exists():
            files = list(dir_path.glob('*'))
            results['datasets_found'].append({
                'dataset': dir_name,
                'files': [f.name for f in files],
                'n_files': len(files)
            })
            results['validation_status'][dir_name] = 'FOUND'
        else:
            results['validation_status'][dir_name] = 'MISSING'
    
    return results


def main():
    """Main function to run comprehensive data validation."""
    print("=" * 80)
    print("COMPREHENSIVE DATA VALIDATION FOR MECHANISM-GWAS-CAUSAL-GRAPHS")
    print("=" * 80)
    print()
    
    # Paths
    data_dir = project_root / 'data'
    processed_dir = data_dir / 'processed'
    raw_dir = data_dir / 'raw'
    manifests_dir = data_dir / 'manifests'
    
    results = {
        'calibration_analysis': {},
        'benchmark_analysis': {},
        'replication_analysis': {},
        'mechanism_graph_analysis': {},
        'locus_summary_analysis': {},
        'manifest_validation': {},
        'overall_validation': {}
    }
    
    # 1. Analyze Calibration Data
    print("1. CALIBRATION ANALYSIS")
    print("-" * 40)
    
    try:
        calibration_metrics = load_tsv(processed_dir / 'calibration' / 'calibration_metrics.tsv')
        reliability_data = load_tsv(processed_dir / 'calibration' / 'reliability_diagram_data.tsv')
        
        calibration_results = analyze_calibration_metrics(calibration_metrics)
        reliability_results = compute_calibration_statistics(reliability_data)
        
        results['calibration_analysis'] = {
            'metrics': calibration_results,
            'reliability': reliability_results
        }
        
        print(f"   Modules analyzed: {len(reliability_results)}")
        for module, stats in reliability_results.items():
            print(f"   - {module}: ECE={stats['ECE']:.4f} ({stats['calibration_status']})")
        
        # Print manuscript claim validation
        if calibration_results.get('manuscript_claims'):
            print("\n   Manuscript Claims Validation:")
            for claim in calibration_results['manuscript_claims']:
                status = "PASS" if claim['validated'] else "FAIL"
                print(f"   - {claim['claim']}: {claim['value']:.3f} [{status}]")
    except Exception as e:
        print(f"   Error in calibration analysis: {e}")
    
    print()
    
    # 2. Analyze Benchmark Data
    print("2. BENCHMARK GENES ANALYSIS")
    print("-" * 40)
    
    try:
        tier1_genes = load_tsv(processed_dir / 'benchmark' / 'tier1_gold_standard_genes.tsv')
        tier2_drugs = load_tsv(processed_dir / 'benchmark' / 'tier2_drug_targets.tsv')
        
        benchmark_results = analyze_benchmark_genes(tier1_genes, tier2_drugs)
        results['benchmark_analysis'] = benchmark_results
        
        t1 = benchmark_results['tier1_analysis']
        t2 = benchmark_results['tier2_analysis']
        
        print(f"   Tier 1 Gold Standard Genes: {t1['total_genes']}")
        print(f"   - Traits covered: {len(t1['traits_covered'])}")
        print(f"   - Anti-leak verified: {t1['anti_leak_verified_pct']}%")
        print(f"   - Training excluded: {t1['training_excluded_pct']}%")
        
        print(f"\n   Tier 2 Drug Targets: {t2['total_drugs']}")
        print(f"   - Unique target genes: {t2['unique_genes']}")
        print(f"   - Approved drugs: {t2['approval_status_distribution'].get('Approved', 0)}")
        
        combined = benchmark_results['combined_analysis']
        print(f"\n   Combined Analysis:")
        print(f"   - Total unique genes: {combined['union']}")
        print(f"   - Overlap T1/T2: {combined['overlap']}")
    except Exception as e:
        print(f"   Error in benchmark analysis: {e}")
    
    print()
    
    # 3. Analyze Replication Data
    print("3. REPLICATION ANALYSIS")
    print("-" * 40)
    
    try:
        replication_data = load_tsv(processed_dir / 'replication' / 'eqtl_catalogue_replication.tsv')
        replication_summary = load_yaml(processed_dir / 'replication' / 'replication_summary.yaml')
        
        replication_results = analyze_replication_data(replication_data, replication_summary)
        results['replication_analysis'] = replication_results
        
        gene_level = replication_results['gene_level']
        print(f"   Genes tested: {gene_level['total_genes']}")
        print(f"   Genes replicated: {gene_level['genes_replicated']}")
        print(f"   Replication rate: {gene_level['replication_rate']}%")
        print(f"   Direction concordance: {gene_level['direction_concordance_rate']}%")
        
        effect_analysis = replication_results.get('effect_size_analysis', {})
        if effect_analysis:
            print(f"\n   Effect Size Analysis:")
            print(f"   - Mean correlation: {effect_analysis.get('mean_correlation', 'N/A')}")
            print(f"   - Median correlation: {effect_analysis.get('median_correlation', 'N/A')}")
        
        claim_validation = replication_results['manuscript_claim_validation']
        print(f"\n   Manuscript Claim Validation:")
        print(f"   - Claimed: 78% replication rate")
        print(f"   - Actual: {claim_validation['actual_replication_rate']}%")
        print(f"   - Status: {claim_validation['status']}")
    except Exception as e:
        print(f"   Error in replication analysis: {e}")
    
    print()
    
    # 4. Analyze Mechanism Graphs
    print("4. MECHANISM GRAPH ANALYSIS")
    print("-" * 40)
    
    try:
        mechanism_graphs = {}
        graphs_dir = processed_dir / 'mechanism_graphs'
        for json_file in graphs_dir.glob('*.json'):
            graph_data = load_json(json_file)
            mechanism_graphs[json_file.stem] = graph_data
        
        graph_results = analyze_mechanism_graphs(mechanism_graphs)
        results['mechanism_graph_analysis'] = graph_results
        
        overall = graph_results['overall']
        print(f"   Loci analyzed: {overall['n_loci']}")
        print(f"   Total mechanism paths: {overall['total_mechanism_paths']}")
        print(f"   Average path probability: {overall['avg_path_probability']:.3f}")
        print(f"   Average colocalization PP.H4: {overall['avg_coloc_pp_h4']:.3f}")
        print(f"   Paths replicated: {overall['paths_replicated']} ({overall['replication_rate']:.1f}%)")
        
        print("\n   Per-Locus Summary:")
        for locus_name, locus_data in graph_results['per_locus'].items():
            print(f"   - {locus_name}:")
            print(f"     Top gene: {locus_data.get('top_gene', 'N/A')} (P={locus_data.get('top_gene_probability', 'N/A')})")
            print(f"     Gold standard: {locus_data.get('is_gold_standard', False)}")
    except Exception as e:
        print(f"   Error in mechanism graph analysis: {e}")
    
    print()
    
    # 5. Analyze Locus Summary
    print("5. LOCUS SUMMARY ANALYSIS")
    print("-" * 40)
    
    try:
        locus_data = load_tsv(processed_dir / 'locus_summary.tsv')
        locus_results = analyze_locus_summary(locus_data)
        results['locus_summary_analysis'] = locus_results
        
        print(f"   Total loci: {locus_results['statistical_summary'].get('total_loci', 0)}")
        
        print("\n   Trait Distribution:")
        for trait, count in locus_results['trait_distribution'].items():
            print(f"   - {trait}: {count}")
        
        print("\n   Evidence Tier Distribution:")
        for tier, count in locus_results['evidence_tier_distribution'].items():
            print(f"   - {tier}: {count}")
        
        print("\n   Validation Status Distribution:")
        for status, count in locus_results['validation_status_distribution'].items():
            print(f"   - {status}: {count}")
        
        stats = locus_results['statistical_summary']
        if 'path_probability' in stats:
            pp = stats['path_probability']
            print(f"\n   Path Probability Statistics:")
            print(f"   - Mean: {pp['mean']}, Median: {pp['median']}")
            print(f"   - Range: [{pp['min']}, {pp['max']}]")
    except Exception as e:
        print(f"   Error in locus summary analysis: {e}")
    
    print()
    
    # 6. Validate Raw Data Manifests
    print("6. RAW DATA MANIFEST VALIDATION")
    print("-" * 40)
    
    try:
        manifest_results = validate_gwas_sumstats_manifest({}, raw_dir)
        results['manifest_validation'] = manifest_results
        
        print("   GWAS Summary Statistics:")
        for dataset in manifest_results['datasets_found']:
            print(f"   - {dataset['dataset']}: {dataset['n_files']} files")
            for f in dataset['files'][:3]:  # Show first 3 files
                print(f"     - {f}")
    except Exception as e:
        print(f"   Error in manifest validation: {e}")
    
    print()
    
    # 7. Overall Validation Summary
    print("=" * 80)
    print("OVERALL VALIDATION SUMMARY")
    print("=" * 80)
    
    validation_passed = []
    validation_failed = []
    
    # Check calibration claims
    if results.get('calibration_analysis', {}).get('metrics', {}).get('manuscript_claims'):
        for claim in results['calibration_analysis']['metrics']['manuscript_claims']:
            if claim['validated']:
                validation_passed.append(claim['claim'])
            else:
                validation_failed.append(claim['claim'])
    
    # Check replication claims
    if results.get('replication_analysis', {}).get('manuscript_claim_validation', {}).get('exceeds_claim'):
        validation_passed.append('Replication rate >= 78%')
    elif 'replication_analysis' in results:
        validation_failed.append('Replication rate >= 78%')
    
    results['overall_validation'] = {
        'claims_validated': validation_passed,
        'claims_failed': validation_failed,
        'total_claims': len(validation_passed) + len(validation_failed),
        'pass_rate': len(validation_passed) / (len(validation_passed) + len(validation_failed)) * 100 if (validation_passed or validation_failed) else 0
    }
    
    print("\n   CLAIMS VALIDATED:")
    for claim in validation_passed:
        print(f"   [PASS] {claim}")
    
    if validation_failed:
        print("\n   CLAIMS NOT VALIDATED:")
        for claim in validation_failed:
            print(f"   [FAIL] {claim}")
    
    print(f"\n   Overall: {len(validation_passed)}/{len(validation_passed) + len(validation_failed)} claims validated")
    print(f"   Pass rate: {results['overall_validation']['pass_rate']:.1f}%")
    
    # Save results
    output_file = project_root / 'comprehensive_validation_results.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    
    print(f"\n   Results saved to: {output_file}")
    print("\n" + "=" * 80)
    print("VALIDATION COMPLETE")
    print("=" * 80)
    
    return results


if __name__ == '__main__':
    main()
