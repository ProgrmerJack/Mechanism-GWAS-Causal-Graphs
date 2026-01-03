#!/usr/bin/env python3
"""
Download Ensembl gene TSS positions using the Ensembl REST API.

For all protein-coding genes in human (GRCh38), retrieves:
- Gene symbol
- Ensembl gene ID
- Chromosome
- TSS position (start for + strand, end for - strand)
- Strand
- Gene biotype

Usage:
    python download_ensembl_tss.py

Output:
    data/external/ensembl/ensembl_tss_grch38.tsv
"""

import requests
import time
import sys
from pathlib import Path
import gzip

def get_all_protein_coding_genes():
    """
    Get all protein-coding genes from Ensembl BioMart.
    
    Uses the Ensembl REST API to retrieve gene information.
    Rate limits: Max 15 requests/second.
    """
    server = "https://rest.ensembl.org"
    
    # Get all genes for human GRCh38
    print("Fetching all protein-coding genes from Ensembl...")
    
    # Use the overlap endpoint to get all genes on each chromosome
    # This is more reliable than trying to get all genes at once
    chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
    
    all_genes = []
    
    for chrom in chromosomes:
        print(f"  Fetching chromosome {chrom}...")
        
        # Use overlap endpoint to get genes on this chromosome
        ext = f"/overlap/region/human/{chrom}:1-300000000?feature=gene;biotype=protein_coding"
        
        try:
            r = requests.get(server + ext, 
                           headers={"Content-Type": "application/json"},
                           timeout=30)
            r.raise_for_status()
            
            genes = r.json()
            print(f"    Found {len(genes)} protein-coding genes")
            
            all_genes.extend(genes)
            
            # Rate limiting
            time.sleep(0.1)
            
        except Exception as e:
            print(f"    Error fetching chromosome {chrom}: {e}")
            continue
    
    print(f"\nTotal protein-coding genes: {len(all_genes)}")
    return all_genes


def get_gene_details_batch(gene_ids, batch_size=200):
    """
    Get detailed gene information for a batch of genes.
    
    Args:
        gene_ids: List of Ensembl gene IDs
        batch_size: Number of genes per batch (max 1000)
    
    Returns:
        Dictionary mapping gene ID to gene details
    """
    server = "https://rest.ensembl.org"
    ext = "/lookup/id"
    
    gene_details = {}
    
    # Process in batches
    for i in range(0, len(gene_ids), batch_size):
        batch = gene_ids[i:i+batch_size]
        
        print(f"  Batch {i//batch_size + 1}/{(len(gene_ids)-1)//batch_size + 1} " +
              f"({len(batch)} genes)...")
        
        try:
            # POST request for batch lookup
            r = requests.post(server + ext,
                            headers={"Content-Type": "application/json",
                                   "Accept": "application/json"},
                            json={"ids": batch},
                            timeout=30)
            r.raise_for_status()
            
            batch_details = r.json()
            gene_details.update(batch_details)
            
            # Rate limiting
            time.sleep(0.1)
            
        except Exception as e:
            print(f"    Error fetching batch: {e}")
            # Try individual lookups for this batch
            for gene_id in batch:
                try:
                    r = requests.get(server + "/lookup/id/" + gene_id,
                                   headers={"Content-Type": "application/json"},
                                   timeout=10)
                    r.raise_for_status()
                    gene_details[gene_id] = r.json()
                    time.sleep(0.1)
                except:
                    continue
    
    return gene_details


def calculate_tss(gene_info):
    """
    Calculate TSS position from gene info.
    
    TSS is:
    - Start position for + strand genes
    - End position for - strand genes
    """
    if gene_info['strand'] == 1:
        return gene_info['start']
    else:
        return gene_info['end']


def main():
    """Main execution function."""
    
    # Setup output directory
    output_dir = Path("C:/Users/Jack0/GitHub/Mechanism-GWAS-Causal-Graphs/data/external/ensembl")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "ensembl_tss_grch38.tsv"
    
    print("="*80)
    print("Ensembl TSS Download Script")
    print("="*80)
    print(f"Output: {output_file}")
    print()
    
    # Get all protein-coding genes
    genes = get_all_protein_coding_genes()
    
    if not genes:
        print("ERROR: No genes retrieved!")
        sys.exit(1)
    
    # Extract TSS information from the overlap data
    print("\nProcessing gene information...")
    
    tss_data = []
    for gene in genes:
        try:
            # Calculate TSS
            tss_pos = calculate_tss(gene)
            
            # Get external_name (gene symbol) if available
            symbol = gene.get('external_name', gene['id'])
            
            tss_data.append({
                'gene_symbol': symbol,
                'ensembl_id': gene['id'],
                'chromosome': gene['seq_region_name'],
                'tss_position': tss_pos,
                'start': gene['start'],
                'end': gene['end'],
                'strand': '+' if gene['strand'] == 1 else '-',
                'biotype': gene.get('biotype', 'protein_coding')
            })
        except Exception as e:
            print(f"  Warning: Could not process {gene.get('id', 'unknown')}: {e}")
            continue
    
    print(f"Processed {len(tss_data)} genes")
    
    # Write to TSV
    print(f"\nWriting to {output_file}...")
    
    with open(output_file, 'w') as f:
        # Header
        f.write("gene_symbol\tensembl_id\tchromosome\ttss_position\tstart\tend\tstrand\tbiotype\n")
        
        # Data rows
        for gene in sorted(tss_data, key=lambda x: (x['chromosome'], x['tss_position'])):
            f.write(f"{gene['gene_symbol']}\t{gene['ensembl_id']}\t" +
                   f"{gene['chromosome']}\t{gene['tss_position']}\t" +
                   f"{gene['start']}\t{gene['end']}\t" +
                   f"{gene['strand']}\t{gene['biotype']}\n")
    
    # Also create gzipped version
    print("Creating gzipped version...")
    with open(output_file, 'rb') as f_in:
        with gzip.open(str(output_file) + '.gz', 'wb') as f_out:
            f_out.write(f_in.read())
    
    print()
    print("="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total genes: {len(tss_data)}")
    print(f"Output file: {output_file}")
    print(f"Compressed: {output_file}.gz")
    
    # Show sample of genes
    print("\nSample genes:")
    for gene in tss_data[:10]:
        print(f"  {gene['gene_symbol']:15s} chr{gene['chromosome']:3s} " +
              f"TSS={gene['tss_position']:11,d} {gene['strand']}")
    
    print("\nDone!")


if __name__ == "__main__":
    main()
