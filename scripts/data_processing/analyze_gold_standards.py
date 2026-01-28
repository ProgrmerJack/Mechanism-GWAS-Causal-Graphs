#!/usr/bin/env python3
"""
Analyze the Open Targets gold standards to extract experimental evidence entries.
These entries are critical for building RegulatoryBench v3 with experimental perturbation evidence.
"""

import json
from pathlib import Path
from collections import Counter
import os

# Change to project root
os.chdir(Path(__file__).parent.parent)

# Load the gold standards - try new path first, then old path
gs_path = Path("data/external/gold_standards/gold_standards/processed/gwas_gold_standards.191108.json")
if not gs_path.exists():
    gs_path = Path("data/external/open_targets/genetics-gold-standards/gold_standards/processed/gwas_gold_standards.191108.json")

with open(gs_path, 'r') as f:
    data = json.load(f)

print(f"Total gold standard entries: {len(data)}")

# Count evidence classes
evidence_classes = Counter()
confidence_levels = Counter()
high_conf_experimental = []
all_experimental = []

for entry in data:
    gs_info = entry.get('gold_standard_info', {})
    evidences = gs_info.get('evidence', [])
    
    for ev in evidences:
        ec = ev.get('class', 'unknown')
        conf = ev.get('confidence', 'unknown')
        
        evidence_classes[ec] += 1
        confidence_levels[conf] += 1
        
        # Collect all experimental entries
        if 'experimental' in ec.lower():
            exp_entry = {
                'gene_id': gs_info.get('gene_id'),
                'rsid': entry.get('sentinel_variant', {}).get('rsid'),
                'chr': entry.get('sentinel_variant', {}).get('locus_GRCh38', {}).get('chromosome'),
                'pos': entry.get('sentinel_variant', {}).get('locus_GRCh38', {}).get('position'),
                'trait': entry.get('trait_info', {}).get('reported_trait_name'),
                'description': ev.get('description', ''),
                'pubmed': ev.get('pubmed_id', ''),
                'confidence': conf,
                'curated_by': ev.get('curated_by', ''),
                'source': ev.get('source', '')
            }
            all_experimental.append(exp_entry)
            
            if conf == 'High':
                high_conf_experimental.append(exp_entry)

print("\n" + "="*60)
print("EVIDENCE CLASS DISTRIBUTION")
print("="*60)
for ec, count in evidence_classes.most_common():
    print(f"  {ec}: {count}")

print("\n" + "="*60)
print("CONFIDENCE LEVEL DISTRIBUTION")
print("="*60)
for conf, count in confidence_levels.most_common():
    print(f"  {conf}: {count}")

print("\n" + "="*60)
print(f"FUNCTIONAL EXPERIMENTAL ENTRIES (Total: {len(all_experimental)}, High-conf: {len(high_conf_experimental)})")
print("="*60)

# Show high-confidence experimental entries
print("\nHigh-confidence functional experimental entries:")
for i, e in enumerate(high_conf_experimental[:30]):
    print(f"\n{i+1}. Gene: {e['gene_id']}, rsID: {e['rsid']}")
    print(f"   Trait: {e['trait']}")
    print(f"   PMID: {e['pubmed']}")
    desc = e['description'][:150] if len(e['description']) > 150 else e['description']
    print(f"   Description: {desc}...")

# Check for CRISPR, reporter assay, MPRA keywords
print("\n" + "="*60)
print("KEYWORD ANALYSIS IN EXPERIMENTAL DESCRIPTIONS")
print("="*60)

keywords = ['CRISPR', 'reporter', 'luciferase', 'MPRA', 'knockdown', 'knockout', 
            'siRNA', 'shRNA', 'overexpression', 'mouse', 'mice', 'zebrafish', 
            'perturbation', 'deletion', 'editing', 'gene editing']

keyword_counts = Counter()
keyword_examples = {}

for e in all_experimental:
    desc = e['description'].lower()
    for kw in keywords:
        if kw.lower() in desc:
            keyword_counts[kw] += 1
            if kw not in keyword_examples:
                keyword_examples[kw] = e

for kw, count in keyword_counts.most_common():
    print(f"  {kw}: {count}")
    ex = keyword_examples.get(kw, {})
    desc = ex.get('description', '')[:100] if ex else ''
    print(f"    Example: {ex.get('gene_id', 'N/A')} - {desc}...")

# Export experimental entries for RegulatoryBench
output_path = Path("data/processed/experimental_gold_standards.json")
output_path.parent.mkdir(parents=True, exist_ok=True)

output_data = {
    'metadata': {
        'source': 'opentargets/genetics-gold-standards',
        'file': 'gwas_gold_standards.191108.json',
        'filter': 'functional experimental evidence only',
        'total_entries': len(all_experimental),
        'high_confidence_entries': len(high_conf_experimental)
    },
    'all_experimental': all_experimental,
    'high_confidence_experimental': high_conf_experimental
}

with open(output_path, 'w') as f:
    json.dump(output_data, f, indent=2)

print(f"\n\nExported experimental entries to: {output_path}")
print(f"Total: {len(all_experimental)}, High-confidence: {len(high_conf_experimental)}")
