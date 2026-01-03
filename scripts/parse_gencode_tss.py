#!/usr/bin/env python3
"""
Parse GENCODE GTF file to extract TSS positions for protein-coding genes.

Creates a TSV file with columns:
- gene_symbol: HGNC gene symbol
- chromosome: Chromosome (1-22, X, Y, MT)
- tss_position: Transcription start site (start for + strand, end for - strand)
- strand: Strand (+/-)
- gene_type: Gene biotype (protein_coding)
"""

import gzip
import sys
from pathlib import Path
from collections import defaultdict

def parse_gtf_attributes(attr_string):
    """Parse GTF attribute field into a dictionary."""
    attrs = {}
    for item in attr_string.strip(';').split(';'):
        item = item.strip()
        if not item:
            continue
        try:
            key, value = item.split(' ', 1)
            attrs[key] = value.strip('"')
        except ValueError:
            # Skip malformed attributes
            continue
    return attrs

def extract_tss_from_gtf(gtf_path, output_path):
    """
    Extract TSS positions for all protein-coding genes from GENCODE GTF.
    
    TSS is defined as:
    - For + strand genes: gene start position
    - For - strand genes: gene end position
    
    Uses the 'gene' feature rows to get gene-level information.
    """
    print("=" * 80)
    print("GENCODE GTF TSS Extraction")
    print("=" * 80)
    print(f"Input:  {gtf_path}")
    print(f"Output: {output_path}")
    print()
    
    # Dictionary to store gene TSS positions
    # Key: gene_symbol, Value: (chromosome, tss_position, strand, gene_type)
    gene_tss = {}
    
    # Track statistics
    total_genes = 0
    protein_coding = 0
    duplicates = 0
    
    print("Parsing GTF file...")
    
    # Open file (handles .gz automatically)
    if str(gtf_path).endswith('.gz'):
        file_handle = gzip.open(gtf_path, 'rt')
    else:
        file_handle = open(gtf_path, 'r')
    
    with file_handle as f:
        for line in f:
            # Skip comments
            if line.startswith('#'):
                continue
            
            # Parse GTF line
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            seqname = fields[0]  # Chromosome
            source = fields[1]
            feature = fields[2]  # Feature type (gene, transcript, exon, etc.)
            start = int(fields[3])
            end = int(fields[4])
            score = fields[5]
            strand = fields[6]
            frame = fields[7]
            attributes = fields[8]
            
            # Only process 'gene' features
            if feature != 'gene':
                continue
            
            total_genes += 1
            
            # Parse attributes
            attrs = parse_gtf_attributes(attributes)
            
            gene_type = attrs.get('gene_type', '')
            gene_name = attrs.get('gene_name', '')
            gene_id = attrs.get('gene_id', '')
            
            # Only keep protein-coding genes
            if gene_type != 'protein_coding':
                continue
            
            protein_coding += 1
            
            # Skip genes without a gene symbol
            if not gene_name:
                continue
            
            # Normalize chromosome names (remove 'chr' prefix if present)
            chrom = seqname.replace('chr', '')
            
            # Calculate TSS
            if strand == '+':
                tss = start
            elif strand == '-':
                tss = end
            else:
                # Unknown strand, skip
                continue
            
            # Check for duplicates
            if gene_name in gene_tss:
                duplicates += 1
                # Keep the first occurrence (usually the canonical transcript)
                continue
            
            # Store TSS position
            gene_tss[gene_name] = (chrom, tss, strand, gene_type)
    
    print(f"  Total genes parsed: {total_genes:,}")
    print(f"  Protein-coding genes: {protein_coding:,}")
    print(f"  Unique gene symbols: {len(gene_tss):,}")
    print(f"  Duplicates skipped: {duplicates:,}")
    print()
    
    # Write output TSV
    print(f"Writing TSS positions to {output_path}...")
    
    with open(output_path, 'w') as out:
        # Write header
        out.write("gene_symbol\tchromosome\ttss_position\tstrand\tgene_type\n")
        
        # Write TSS positions, sorted by gene symbol
        for gene_symbol in sorted(gene_tss.keys()):
            chrom, tss, strand, gene_type = gene_tss[gene_symbol]
            out.write(f"{gene_symbol}\t{chrom}\t{tss}\t{strand}\t{gene_type}\n")
    
    print(f"  Wrote {len(gene_tss):,} protein-coding genes")
    print()
    print("SUCCESS: TSS extraction complete!")
    print("=" * 80)
    
    return gene_tss

def main():
    """Main entry point."""
    # Define paths
    base_dir = Path(__file__).parent.parent
    gtf_path = base_dir / 'data' / 'external' / 'ensembl' / 'gencode.v38.annotation.gtf.gz'
    output_path = base_dir / 'data' / 'external' / 'ensembl' / 'gencode_tss_grch38.tsv'
    
    # Check if GTF exists
    if not gtf_path.exists():
        print(f"ERROR: GTF file not found: {gtf_path}")
        print("Please download GENCODE GTF first:")
        print("  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz")
        sys.exit(1)
    
    # Extract TSS positions
    extract_tss_from_gtf(gtf_path, output_path)

if __name__ == '__main__':
    main()
