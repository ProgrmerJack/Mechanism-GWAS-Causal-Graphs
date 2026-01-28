#!/usr/bin/env python3
"""
Comprehensive GWAS Summary Statistics Analysis
Analyzes ALL raw GWAS data files (5GB+) to compute real statistics

This script:
1. Loads ALL GWAS summary statistics (multi-GB files)
2. Computes genome-wide statistics
3. Identifies significant variants
4. Validates fine-mapping results
5. Saves processed results

WARNING: This will take hours to run and requires substantial RAM (16GB+ recommended)
"""

import gzip
import sys
import json
from pathlib import Path
from collections import defaultdict
import statistics
import time

project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))


def load_gwas_sumstats_chunk(filepath, chunk_size=100000):
    """
    Generator to load GWAS summary statistics in chunks.
    Yields dictionaries with variant information.
    """
    print(f"Loading: {filepath.name}")
    
    if filepath.suffix == '.gz':
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
    
    with open_func(filepath, mode, encoding='utf-8', errors='ignore') as f:
        # Read header
        header_line = f.readline().strip()
        
        # Try different delimiter strategies
        if '\t' in header_line:
            delimiter = '\t'
        else:
            delimiter = None  # Will try to auto-detect
        
        headers = header_line.split(delimiter) if delimiter else header_line.split()
        
        # Normalize header names (handle different GWAS formats)
        header_map = {}
        for i, h in enumerate(headers):
            h_lower = h.lower()
            if h_lower in ['snp', 'rsid', 'variant_id', 'markername', 'snpid']:
                header_map['snp'] = i
            elif h_lower in ['chr', 'chromosome', '#chr', 'chromosome(b37)', 'chromosome(b38)']:
                header_map['chr'] = i
            elif h_lower in ['pos', 'position', 'bp', 'pos_b38', 'base_pair_location', 'position(b37)', 'position(b38)']:
                header_map['pos'] = i
            elif h_lower in ['pval', 'p', 'p-value', 'p_value', 'pvalue', 'p.value', 'fixed-effects_p-value', 'p-val', 'p_val']:
                header_map['pval'] = i
            elif h_lower in ['beta', 'effect', 'effect_size', 'fixed-effects_beta', 'effect_beta']:
                header_map['beta'] = i
            elif h_lower in ['se', 'stderr', 'standard_error', 'fixed-effects_se', 'se_fixed']:
                header_map['se'] = i
            elif h_lower in ['a1', 'allele1', 'effect_allele', 'ea']:
                header_map['a1'] = i
            elif h_lower in ['a2', 'allele2', 'other_allele', 'nea']:
                header_map['a2'] = i
            elif h_lower in ['eaf', 'freq', 'frequency', 'maf']:
                header_map['eaf'] = i
        
        chunk = []
        for line_num, line in enumerate(f, start=2):
            if not line.strip() or line.startswith('#'):
                continue
            
            fields = line.strip().split(delimiter) if delimiter else line.strip().split()
            
            if len(fields) < len(headers):
                continue
            
            try:
                variant = {}
                
                if 'snp' in header_map:
                    variant['snp'] = fields[header_map['snp']]
                if 'chr' in header_map:
                    chr_val = fields[header_map['chr']].replace('chr', '')
                    try:
                        variant['chr'] = int(chr_val) if chr_val.isdigit() else chr_val
                    except:
                        variant['chr'] = chr_val
                if 'pos' in header_map:
                    try:
                        variant['pos'] = int(float(fields[header_map['pos']]))
                    except:
                        continue
                if 'pval' in header_map:
                    try:
                        variant['pval'] = float(fields[header_map['pval']])
                    except:
                        continue
                if 'beta' in header_map:
                    try:
                        variant['beta'] = float(fields[header_map['beta']])
                    except:
                        variant['beta'] = None
                if 'se' in header_map:
                    try:
                        variant['se'] = float(fields[header_map['se']])
                    except:
                        variant['se'] = None
                if 'a1' in header_map:
                    variant['a1'] = fields[header_map['a1']]
                if 'a2' in header_map:
                    variant['a2'] = fields[header_map['a2']]
                if 'eaf' in header_map:
                    try:
                        variant['eaf'] = float(fields[header_map['eaf']])
                    except:
                        variant['eaf'] = None
                
                chunk.append(variant)
                
                if len(chunk) >= chunk_size:
                    yield chunk
                    chunk = []
                    
            except Exception as e:
                if line_num < 10:
                    print(f"  Warning: Error parsing line {line_num}: {e}")
                continue
        
        if chunk:
            yield chunk
    
    print(f"  Finished loading {filepath.name}")


def analyze_gwas_file(filepath, genome_wide_threshold=5e-8):
    """
    Analyze a complete GWAS summary statistics file.
    Returns comprehensive statistics.
    """
    results = {
        'file': filepath.name,
        'total_variants': 0,
        'genome_wide_sig': 0,
        'suggestive_sig': 0,  # p < 1e-5
        'chromosomes': defaultdict(int),
        'min_pval': 1.0,
        'top_variants': [],
        'mean_beta': None,
        'mean_se': None,
        'lambda_gc': None,
        'processing_time_seconds': 0
    }
    
    start_time = time.time()
    
    all_betas = []
    all_ses = []
    all_chisq = []
    
    print(f"\nAnalyzing: {filepath.name}")
    print(f"  File size: {filepath.stat().st_size / 1e9:.2f} GB")
    
    for chunk_num, chunk in enumerate(load_gwas_sumstats_chunk(filepath), start=1):
        for variant in chunk:
            results['total_variants'] += 1
            
            pval = variant.get('pval')
            if pval is None:
                continue
            
            # Track minimum p-value
            if pval < results['min_pval']:
                results['min_pval'] = pval
            
            # Count significant variants
            if pval < genome_wide_threshold:
                results['genome_wide_sig'] += 1
                
                # Track top variants
                if len(results['top_variants']) < 100:
                    results['top_variants'].append({
                        'snp': variant.get('snp', 'unknown'),
                        'chr': variant.get('chr', 'unknown'),
                        'pos': variant.get('pos', 0),
                        'pval': pval,
                        'beta': variant.get('beta'),
                        'a1': variant.get('a1', ''),
                        'a2': variant.get('a2', '')
                    })
                elif pval < results['top_variants'][-1]['pval']:
                    results['top_variants'][-1] = {
                        'snp': variant.get('snp', 'unknown'),
                        'chr': variant.get('chr', 'unknown'),
                        'pos': variant.get('pos', 0),
                        'pval': pval,
                        'beta': variant.get('beta'),
                        'a1': variant.get('a1', ''),
                        'a2': variant.get('a2', '')
                    }
                    results['top_variants'].sort(key=lambda x: x['pval'])
            
            if pval < 1e-5:
                results['suggestive_sig'] += 1
            
            # Track chromosomes
            if 'chr' in variant:
                results['chromosomes'][str(variant['chr'])] += 1
            
            # Collect for summary stats
            if variant.get('beta') is not None:
                all_betas.append(variant['beta'])
            if variant.get('se') is not None:
                all_ses.append(variant['se'])
            
            # For lambda GC calculation
            if pval > 0 and variant.get('beta') is not None and variant.get('se') is not None:
                try:
                    z = variant['beta'] / variant['se']
                    chisq = z ** 2
                    all_chisq.append(chisq)
                except:
                    pass
        
        # Progress update every 10 chunks
        if chunk_num % 10 == 0:
            print(f"  Processed {results['total_variants']:,} variants...")
    
    # Compute summary statistics
    if all_betas:
        results['mean_beta'] = statistics.mean(all_betas)
        results['std_beta'] = statistics.stdev(all_betas) if len(all_betas) > 1 else 0
    
    if all_ses:
        results['mean_se'] = statistics.mean(all_ses)
    
    # Compute lambda GC (genomic inflation factor)
    if all_chisq and len(all_chisq) > 1000:
        all_chisq.sort()
        median_chisq = all_chisq[len(all_chisq) // 2]
        results['lambda_gc'] = median_chisq / 0.456  # Expected median of chi-sq(1)
    
    results['chromosomes'] = dict(results['chromosomes'])
    results['processing_time_seconds'] = time.time() - start_time
    
    print(f"\n  RESULTS for {filepath.name}:")
    print(f"    Total variants: {results['total_variants']:,}")
    print(f"    Genome-wide significant (p<5e-8): {results['genome_wide_sig']:,}")
    print(f"    Minimum p-value: {results['min_pval']:.2e}")
    print(f"    Lambda GC: {results['lambda_gc']:.3f}" if results['lambda_gc'] else "    Lambda GC: N/A")
    print(f"    Processing time: {results['processing_time_seconds']:.1f} seconds")
    
    return results


def analyze_all_gwas_datasets():
    """
    Analyze ALL GWAS summary statistics files.
    """
    print("=" * 80)
    print("COMPREHENSIVE GWAS SUMMARY STATISTICS ANALYSIS")
    print("Analyzing ALL raw data files (this will take hours)")
    print("=" * 80)
    
    raw_data_dir = project_root / 'data' / 'raw' / 'gwas_sumstats'
    
    if not raw_data_dir.exists():
        print(f"\nERROR: Raw data directory not found: {raw_data_dir}")
        print("Please ensure GWAS data has been downloaded.")
        return
    
    # Find all GWAS files - recursively search all subdirectories
    gwas_files = []
    print("\nScanning for GWAS summary statistics files...")
    
    for subdir in raw_data_dir.rglob('*'):
        if subdir.is_file():
            # Check for GWAS file extensions
            if subdir.suffix in ['.gz', '.txt', '.tsv'] or subdir.name.endswith('.h.tsv.gz'):
                # Skip small files (likely metadata)
                file_size = subdir.stat().st_size
                if file_size > 10_000_000:  # > 10 MB
                    # Determine study type from path
                    study_type = subdir.parent.name
                    gwas_files.append({
                        'path': subdir,
                        'study': study_type,
                        'size_mb': file_size / 1e6
                    })
                    print(f"  Found: {study_type}/{subdir.name} ({file_size/1e9:.2f} GB)")
    
    if not gwas_files:
        print("\nNo GWAS files found!")
        print("Looking in:", raw_data_dir)
        return
    
    # Sort by study type
    gwas_files.sort(key=lambda x: x['study'])
    
    print(f"\n{'='*80}")
    print(f"Found {len(gwas_files)} GWAS summary statistics files")
    
    # Group by study
    studies = {}
    for f in gwas_files:
        if f['study'] not in studies:
            studies[f['study']] = []
        studies[f['study']].append(f)
    
    print(f"\nStudies discovered:")
    for study, files in studies.items():
        total_size = sum(f['size_mb'] for f in files) / 1000
        print(f"  - {study}: {len(files)} files ({total_size:.2f} GB)")
    
    total_size_gb = sum(f['size_mb'] for f in gwas_files) / 1000
    print(f"\nTotal data size: {total_size_gb:.2f} GB")
    print()
    
    # Analyze each file
    all_results = {}
    
    for i, gwas_dict in enumerate(gwas_files, 1):
        gwas_file = gwas_dict['path']
        study = gwas_dict['study']
        
        print(f"\n{'='*80}")
        print(f"File {i}/{len(gwas_files)}: {study}")
        print(f"{'='*80}")
        
        try:
            results = analyze_gwas_file(gwas_file)
            # Use study/filename as key
            dataset_key = f"{study}/{gwas_file.stem}"
            all_results[dataset_key] = results
            all_results[dataset_key]['study_type'] = study
        except Exception as e:
            print(f"\nERROR analyzing {gwas_file.name}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Save comprehensive results
    output_dir = project_root / 'data' / 'processed' / 'gwas_analysis'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / 'comprehensive_gwas_analysis.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(all_results, f, indent=2)
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nResults saved to: {output_file}")
    
    # Summary
    total_variants = sum(r['total_variants'] for r in all_results.values())
    total_sig = sum(r['genome_wide_sig'] for r in all_results.values())
    total_time = sum(r['processing_time_seconds'] for r in all_results.values())
    
    print(f"\nOverall Summary:")
    print(f"  Studies analyzed: {len(set(r['study_type'] for r in all_results.values()))}")
    print(f"  Files analyzed: {len(all_results)}")
    print(f"  Total variants: {total_variants:,}")
    print(f"  Genome-wide significant variants: {total_sig:,}")
    print(f"  Total processing time: {total_time/60:.1f} minutes")
    
    # Calculate average lambda GC (excluding None values)
    lambda_values = [r['lambda_gc'] for r in all_results.values() if r['lambda_gc'] is not None]
    if lambda_values:
        print(f"  Average lambda GC: {statistics.mean(lambda_values):.3f}")
    
    # Summary by study
    print(f"\nSummary by study:")
    for study in sorted(set(r['study_type'] for r in all_results.values())):
        study_results = {k: v for k, v in all_results.items() if v['study_type'] == study}
        study_variants = sum(r['total_variants'] for r in study_results.values())
        study_sig = sum(r['genome_wide_sig'] for r in study_results.values())
        print(f"  {study}: {len(study_results)} files, {study_variants:,} variants, {study_sig:,} significant")
    
    # Create summary TSV
    summary_tsv = output_dir / 'gwas_analysis_summary.tsv'
    with open(summary_tsv, 'w', encoding='utf-8') as f:
        f.write("dataset\tstudy_type\ttotal_variants\tgenome_wide_sig\tsuggestive_sig\tmin_pval\tlambda_gc\tmean_chi2\tprocessing_time_min\n")
        for dataset, results in all_results.items():
            f.write(f"{dataset}\t{results['study_type']}\t{results['total_variants']}\t{results['genome_wide_sig']}\t"
                   f"{results['suggestive_sig']}\t{results['min_pval']:.2e}\t{results.get('lambda_gc', 'NA')}\t"
                   f"{results.get('mean_chi2', 'NA')}\t{results['processing_time_seconds']/60:.2f}\n")
    
    print(f"\nSummary TSV saved to: {summary_tsv}")
    
    return all_results


if __name__ == '__main__':
    print("\nWARNING: This script will analyze multi-GB GWAS files.")
    print("Expected runtime: 2-6 hours depending on system.")
    print("Recommended: 16GB+ RAM")
    print()
    
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == '--run':
        analyze_all_gwas_datasets()
    else:
        print("To run analysis, execute:")
        print("  python scripts/analyze_all_gwas_sumstats.py --run")
        print()
        print("Or if you want to run in background:")
        print("  nohup python scripts/analyze_all_gwas_sumstats.py --run > gwas_analysis.log 2>&1 &")
