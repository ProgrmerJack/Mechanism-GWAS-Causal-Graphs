#!/usr/bin/env python3
"""
Figure Quality Audit Script
Analyzes all manuscript figures for Nature Genetics compliance.
"""

import fitz  # PyMuPDF
from pathlib import Path
import json

PROJECT_ROOT = Path(__file__).parent.parent
FIGURES_DIR = PROJECT_ROOT / "manuscript" / "figures"

MAIN_FIGURES = [
    'fig1_overview.pdf',
    'fig2_bridge.pdf', 
    'fig3_benchmark.pdf',
    'fig4_calibration.pdf',
    'fig5_pqtl_validation.pdf',
    'fig6_examples.pdf',
]

ED_FIGURES = [
    'ed_fig1_datasets.pdf',
    'ed_fig2_multicausal.pdf',
    'ed_fig3_benchmark_gene_provenance.pdf',
    'ed_fig4_eqtl_catalogue_tissue_matching.pdf',
    'ed_fig5_correlation_correction_validation.pdf',
    'ed_fig6_out-of-domain_performance_details.pdf',
    'ed_fig7_failure_mode_examples.pdf',
    'ed_fig8_bootstrap_confidence_intervals.pdf',
    'ed_fig9_negative_controls.pdf',
    'ed_fig10_reliability_diagrams.pdf',
]

# Nature Genetics requirements
NG_REQUIREMENTS = {
    'single_column_width': 3.5,  # inches (89mm)
    'double_column_width': 7.2,  # inches (183mm)
    'max_height': 9.0,  # inches
    'min_font_size': 5,  # points
    'recommended_font_size': 7,  # points
    'max_font_size': 8,  # points for labels
    'line_width_min': 0.5,
    'dpi_min': 300,
}

def analyze_pdf(filepath):
    """Analyze a PDF figure and return quality metrics."""
    try:
        doc = fitz.open(str(filepath))
        page = doc[0]
        rect = page.rect
        
        # Dimensions
        width_inches = rect.width / 72
        height_inches = rect.height / 72
        
        # Text extraction
        text = page.get_text()
        text_blocks = page.get_text("blocks")
        
        # Image analysis
        images = page.get_images()
        
        # Font analysis  
        fonts = set()
        text_dict = page.get_text("dict")
        for block in text_dict.get("blocks", []):
            if "lines" in block:
                for line in block["lines"]:
                    for span in line["spans"]:
                        fonts.add((span.get("font", "unknown"), span.get("size", 0)))
        
        doc.close()
        
        return {
            'exists': True,
            'width_inches': width_inches,
            'height_inches': height_inches,
            'has_text': len(text.strip()) > 0,
            'text_length': len(text),
            'num_text_blocks': len(text_blocks),
            'num_images': len(images),
            'fonts': list(fonts)[:10],  # First 10 fonts
            'text_sample': text[:200].replace('\n', ' ').strip(),
        }
    except Exception as e:
        return {'exists': False, 'error': str(e)}


def grade_figure(metrics):
    """Grade figure quality based on Nature Genetics requirements."""
    issues = []
    score = 100
    
    if not metrics.get('exists'):
        return 0, ['FILE NOT FOUND']
    
    width = metrics['width_inches']
    height = metrics['height_inches']
    
    # Size check
    if width > NG_REQUIREMENTS['double_column_width'] + 0.5:
        issues.append(f"TOO WIDE: {width:.1f}\" (max {NG_REQUIREMENTS['double_column_width']}\")")
        score -= 20
    
    if height > NG_REQUIREMENTS['max_height']:
        issues.append(f"TOO TALL: {height:.1f}\" (max {NG_REQUIREMENTS['max_height']}\")")
        score -= 15
    
    # Font check
    fonts = metrics.get('fonts', [])
    for font_name, font_size in fonts:
        if font_size > 0 and font_size < NG_REQUIREMENTS['min_font_size']:
            issues.append(f"FONT TOO SMALL: {font_size:.1f}pt (min {NG_REQUIREMENTS['min_font_size']}pt)")
            score -= 10
            break
    
    # Check for raster images (bad for publication)
    if metrics.get('num_images', 0) > 0:
        issues.append(f"CONTAINS RASTER: {metrics['num_images']} embedded image(s) - should be vector")
        score -= 15
    
    # Text check
    if not metrics.get('has_text'):
        issues.append("NO EDITABLE TEXT - fonts may be outlined/rasterized")
        score -= 25
    
    # Professional appearance
    if metrics.get('text_length', 0) < 50:
        issues.append("MINIMAL TEXT - may appear unprofessional")
        score -= 10
    
    return max(0, score), issues


def main():
    print("=" * 80)
    print("NATURE GENETICS FIGURE QUALITY AUDIT")
    print("=" * 80)
    
    all_results = {}
    
    print("\n" + "=" * 40)
    print("MAIN FIGURES")
    print("=" * 40)
    
    for fig_name in MAIN_FIGURES:
        fig_path = FIGURES_DIR / fig_name
        metrics = analyze_pdf(fig_path)
        score, issues = grade_figure(metrics)
        all_results[fig_name] = {'metrics': metrics, 'score': score, 'issues': issues}
        
        print(f"\n📊 {fig_name}")
        print(f"   Score: {score}/100 {'✅' if score >= 80 else '⚠️' if score >= 60 else '❌'}")
        
        if metrics.get('exists'):
            print(f"   Size: {metrics['width_inches']:.1f}\" × {metrics['height_inches']:.1f}\"")
            print(f"   Text blocks: {metrics['num_text_blocks']}, Images: {metrics['num_images']}")
            if metrics.get('fonts'):
                font_sizes = [f[1] for f in metrics['fonts'] if f[1] > 0]
                if font_sizes:
                    print(f"   Font sizes: {min(font_sizes):.1f} - {max(font_sizes):.1f} pt")
        
        for issue in issues:
            print(f"   ❗ {issue}")
    
    print("\n" + "=" * 40)
    print("EXTENDED DATA FIGURES")
    print("=" * 40)
    
    for fig_name in ED_FIGURES:
        fig_path = FIGURES_DIR / fig_name
        metrics = analyze_pdf(fig_path)
        score, issues = grade_figure(metrics)
        all_results[fig_name] = {'metrics': metrics, 'score': score, 'issues': issues}
        
        print(f"\n📊 {fig_name}")
        print(f"   Score: {score}/100 {'✅' if score >= 80 else '⚠️' if score >= 60 else '❌'}")
        
        if metrics.get('exists'):
            print(f"   Size: {metrics['width_inches']:.1f}\" × {metrics['height_inches']:.1f}\"")
        
        for issue in issues[:3]:  # Limit to 3 issues for brevity
            print(f"   ❗ {issue}")
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    scores = [r['score'] for r in all_results.values()]
    avg_score = sum(scores) / len(scores) if scores else 0
    
    critical_issues = sum(1 for r in all_results.values() if r['score'] < 60)
    warnings = sum(1 for r in all_results.values() if 60 <= r['score'] < 80)
    good = sum(1 for r in all_results.values() if r['score'] >= 80)
    
    print(f"\nOverall Quality Score: {avg_score:.0f}/100")
    print(f"✅ Good (≥80): {good}")
    print(f"⚠️  Warnings (60-79): {warnings}")
    print(f"❌ Critical (<60): {critical_issues}")
    
    # Most common issues
    all_issues = []
    for r in all_results.values():
        all_issues.extend(r['issues'])
    
    print("\n🔧 Most Common Issues:")
    issue_types = {}
    for issue in all_issues:
        issue_type = issue.split(':')[0]
        issue_types[issue_type] = issue_types.get(issue_type, 0) + 1
    
    for issue_type, count in sorted(issue_types.items(), key=lambda x: -x[1])[:5]:
        print(f"   - {issue_type}: {count} figures")
    
    return all_results


if __name__ == '__main__':
    main()
