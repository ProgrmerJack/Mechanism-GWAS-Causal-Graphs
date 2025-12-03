# Data Manifests for Reproducibility

This directory contains versioned data manifests with cryptographic hashes
for complete reproducibility of all analyses in this project.

## Manifest Files

| File | Description |
|------|-------------|
| `gwas_sumstats.yaml` | GWAS summary statistics sources with MD5/SHA256 hashes |
| `reference_data.yaml` | Reference genome and annotation files |
| `functional_data.yaml` | ENCODE cCREs, ABC scores, PCHi-C contacts |
| `eqtl_data.yaml` | GTEx v8 and eQTL Catalogue data |
| `benchmark_genes.yaml` | Gold-standard gene sets with provenance |

## Verification

To verify data integrity:

```bash
python scripts/verify_manifests.py --manifest data/manifests/
```

## Provenance Tracking

Each manifest includes:
- **Source URL**: Original download location
- **Version**: Data version/release date
- **MD5/SHA256**: Cryptographic hash for integrity
- **Download date**: When data was acquired
- **License**: Data usage terms

## Updates

When updating any data file:
1. Download new version
2. Compute hashes: `sha256sum filename`
3. Update relevant manifest
4. Commit manifest change with version note
