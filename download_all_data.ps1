# Comprehensive Data Download Script for Mechanism-GWAS-Causal-Graphs
# Based on your actual project manifest files (data/manifests/*.yaml)
# Run this script to download all required data for the analysis

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Mechanism-GWAS-Causal-Graphs Data Download" -ForegroundColor Cyan
Write-Host "Based on YOUR project manifest files" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

# Create data directories
$dataDir = "data/raw"
New-Item -ItemType Directory -Force -Path "$dataDir/gwas_catalog" | Out-Null
New-Item -ItemType Directory -Force -Path "$dataDir/encode_ccre" | Out-Null
New-Item -ItemType Directory -Force -Path "$dataDir/gtex_v8" | Out-Null
New-Item -ItemType Directory -Force -Path "$dataDir/eqtl_catalogue" | Out-Null
New-Item -ItemType Directory -Force -Path "$dataDir/ld_reference" | Out-Null
New-Item -ItemType Directory -Force -Path "$dataDir/gene_annotations" | Out-Null
New-Item -ItemType Directory -Force -Path "data/external" | Out-Null

Write-Host "Created directory structure" -ForegroundColor Green
Write-Host ""

# Function to download with progress
function Download-FileWithProgress {
    param(
        [string]$Url,
        [string]$OutFile,
        [string]$Description
    )
    
    Write-Host "Downloading: $Description" -ForegroundColor Yellow
    Write-Host "URL: $Url" -ForegroundColor Gray
    Write-Host "Output: $OutFile" -ForegroundColor Gray
    
    try {
        $ProgressPreference = 'SilentlyContinue'
        Invoke-WebRequest -Uri $Url -OutFile $OutFile -TimeoutSec 300
        Write-Host "✓ Downloaded successfully" -ForegroundColor Green
        Write-Host ""
        return $true
    }
    catch {
        Write-Host "✗ Failed: $_" -ForegroundColor Red
        Write-Host "Please download manually from: $Url" -ForegroundColor Yellow
        Write-Host ""
        return $false
    }
}

# ============================================================================
# 1. ENCODE cCRE Data (Small files, should work)
# ============================================================================
Write-Host "1. ENCODE cCRE Annotations" -ForegroundColor Cyan
Write-Host "----------------------------------------" -ForegroundColor Cyan

$ccreFiles = @{
    "GRCh38-cCREs.bed.gz" = "https://api.wenglab.org/screen/downloads/Seven-Group/V3/GRCh38-cCREs.bed.gz"
    "GRCh38-PLS.bed.gz" = "https://api.wenglab.org/screen/downloads/cCREs/GRCh38-PLS.bed"
    "GRCh38-pELS.bed.gz" = "https://api.wenglab.org/screen/downloads/cCREs/GRCh38-pELS.bed"
    "GRCh38-dELS.bed.gz" = "https://api.wenglab.org/screen/downloads/cCREs/GRCh38-dELS.bed"
    "GRCh38-CTCF.bed.gz" = "https://api.wenglab.org/screen/downloads/cCREs/GRCh38-CTCF-only.bed"
}

foreach ($file in $ccreFiles.GetEnumerator()) {
    $outPath = Join-Path "$dataDir/encode_ccre" $file.Key
    if (-not (Test-Path $outPath)) {
        Download-FileWithProgress -Url $file.Value -OutFile $outPath -Description "cCRE: $($file.Key)"
    } else {
        Write-Host "✓ Already exists: $($file.Key)" -ForegroundColor Green
    }
}

# ============================================================================
# 2. GENCODE Gene Annotations v26 (matches GTEx v8)
# ============================================================================
Write-Host ""
Write-Host "2. GENCODE Gene Annotations v26" -ForegroundColor Cyan
Write-Host "----------------------------------------" -ForegroundColor Cyan
Write-Host "NOTE: Your manifest specifies v26 to match GTEx v8" -ForegroundColor Yellow

$gencodeFiles = @{
    "gencode.v26.annotation.gtf.gz" = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz"
}

foreach ($file in $gencodeFiles.GetEnumerator()) {
    $outPath = Join-Path "$dataDir/gene_annotations" $file.Key
    if (-not (Test-Path $outPath)) {
        Download-FileWithProgress -Url $file.Value -OutFile $outPath -Description $file.Key
    } else {
        Write-Host "✓ Already exists: $($file.Key)" -ForegroundColor Green
    }
}

# ============================================================================
# 3. GTEx v8 eQTL Data (LARGE FILES - may timeout)
# ============================================================================
Write-Host ""
Write-Host "3. GTEx v8 eQTL Data" -ForegroundColor Cyan
Write-Host "----------------------------------------" -ForegroundColor Cyan
Write-Host "Note: These are LARGE files (100MB-500MB each). Downloads may timeout." -ForegroundColor Yellow
Write-Host "If download fails, please download manually from:" -ForegroundColor Yellow
Write-Host "https://gtexportal.org/home/datasets" -ForegroundColor Yellow
Write-Host ""

$tissues = @(
    "Liver",
    "Adipose_Subcutaneous",
    "Adipose_Visceral_Omentum",
    "Heart_Left_Ventricle",
    "Heart_Atrial_Appendage",
    "Artery_Aorta",
    "Artery_Coronary",
    "Artery_Tibial",
    "Pancreas",
    "Muscle_Skeletal",
    "Whole_Blood"
)

$gtexBaseUrl = "https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL"

foreach ($tissue in $tissues) {
    $fileName = "$tissue.v8.signif_variant_gene_pairs.txt.gz"
    $url = "$gtexBaseUrl/$fileName"
    $outPath = Join-Path "$dataDir/gtex_v8" $fileName
    
    if (-not (Test-Path $outPath)) {
        Write-Host "Attempting: $fileName" -ForegroundColor Yellow
        $success = Download-FileWithProgress -Url $url -OutFile $outPath -Description "GTEx eQTL: $tissue"
        if (-not $success) {
            Write-Host "MANUAL DOWNLOAD REQUIRED:" -ForegroundColor Red
            Write-Host "wget $url -O $outPath" -ForegroundColor White
            Write-Host ""
        }
    } else {
        Write-Host "✓ Already exists: $fileName" -ForegroundColor Green
    }
}

# ============================================================================
# 4. GWAS Summary Statistics (VERY LARGE - manual download recommended)
# ============================================================================
Write-Host ""
Write-Host "4. GWAS Catalog Summary Statistics" -ForegroundColor Cyan
Write-Host "----------------------------------------" -ForegroundColor Cyan
Write-Host "These files are VERY LARGE (1-10GB each)." -ForegroundColor Yellow
Write-Host "Automated download not recommended. Please download manually:" -ForegroundColor Yellow
Write-Host ""

$gwasDownloads = @"
# GLGC Lipids (2021)
http://csg.sph.umich.edu/willer/public/glgc-lipids2021/GLGC_EUR_LDL.txt.gz
http://csg.sph.umich.edu/willer/public/glgc-lipids2021/GLGC_EUR_HDL.txt.gz
http://csg.sph.umich.edu/willer/public/glgc-lipids2021/GLGC_EUR_TG.txt.gz
http://csg.sph.umich.edu/willer/public/glgc-lipids2021/GLGC_EUR_TC.txt.gz

# CAD (CARDIoGRAMplusC4D 2022)
http://www.cardiogramplusc4d.org/data-downloads/

# Type 2 Diabetes (DIAGRAM)
https://diagram-consortium.org/downloads.html

# Blood Pressure (ICBP)
https://www.ebi.ac.uk/gwas/studies/GCST006624
https://www.ebi.ac.uk/gwas/studies/GCST006629
"@

Write-Host $gwasDownloads -ForegroundColor White
Write-Host ""

# Save URLs to file for reference
$gwasDownloads | Out-File -FilePath "$dataDir/gwas_catalog/DOWNLOAD_URLS.txt" -Encoding UTF8
Write-Host "✓ Saved download URLs to: $dataDir/gwas_catalog/DOWNLOAD_URLS.txt" -ForegroundColor Green

# ============================================================================
# 5. 1000 Genomes LD Reference (EXTREMELY LARGE - manual download required)
# ============================================================================
Write-Host ""
Write-Host "5. 1000 Genomes Phase 3 LD Reference" -ForegroundColor Cyan
Write-Host "----------------------------------------" -ForegroundColor Cyan
Write-Host "EXTREMELY LARGE FILES (~100GB total)." -ForegroundColor Red
Write-Host "Manual download required using wget or specialized tools." -ForegroundColor Yellow
Write-Host ""

$ldScript = @"
# Download 1000 Genomes Phase 3 VCF files
# European population for LD calculation

BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
OUTPUT_DIR="data/raw/ld_reference"

mkdir -p `$OUTPUT_DIR

# Download chromosome-by-chromosome
for chr in {1..22}; do
    echo "Downloading chromosome `$chr..."
    wget `$BASE_URL/ALL.chr`${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
         -O `$OUTPUT_DIR/chr`${chr}.vcf.gz
    
    wget `$BASE_URL/ALL.chr`${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi \
         -O `$OUTPUT_DIR/chr`${chr}.vcf.gz.tbi
done

# Download sample information
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index \
     -O `$OUTPUT_DIR/sample_info.txt
"@

$ldScript | Out-File -FilePath "$dataDir/ld_reference/download_1000g.sh" -Encoding UTF8
Write-Host "✓ Created download script: $dataDir/ld_reference/download_1000g.sh" -ForegroundColor Green
Write-Host "  Run this on a Linux/Mac system with wget installed" -ForegroundColor Gray

# ============================================================================
# 6. Zenodo Repository
# ============================================================================
Write-Host ""
Write-Host "6. Zenodo Repository (Processed Data)" -ForegroundColor Cyan
Write-Host "----------------------------------------" -ForegroundColor Cyan
Write-Host "DOI: 10.5281/zenodo.17799751" -ForegroundColor White
Write-Host "URL: https://doi.org/10.5281/zenodo.17799751" -ForegroundColor White
Write-Host ""
Write-Host "This contains processed data and may be easier than raw downloads." -ForegroundColor Yellow
Write-Host "Visit the URL above to download the complete dataset." -ForegroundColor Yellow

# ============================================================================
# Summary
# ============================================================================
Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Download Summary" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Small files (cCREs, gene annotations): Attempted automatic download" -ForegroundColor Green
Write-Host "Medium files (GTEx eQTLs): Attempted automatic download (may have failed)" -ForegroundColor Yellow
Write-Host "Large files (GWAS, 1000G): Manual download required - see URLs above" -ForegroundColor Red
Write-Host ""
Write-Host "Next steps:" -ForegroundColor Cyan
Write-Host "1. Check which downloads succeeded" -ForegroundColor White
Write-Host "2. Manually download large files using provided URLs" -ForegroundColor White
Write-Host "3. Or use your Python script: python scripts/download_data.py --config config/config.yaml" -ForegroundColor White
Write-Host ""
Write-Host "Full guide with correct sources: ACTUAL_DATA_SOURCES_GUIDE.md" -ForegroundColor Yellow
Write-Host "Your project repo: https://github.com/ProgrmerJack/Mechanism-GWAS-Causal-Graphs" -ForegroundColor White
