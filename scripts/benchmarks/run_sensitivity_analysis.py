
import sys
from pathlib import Path
import logging
import pandas as pd

# Add project root to path
project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root))

from scripts.run_cs2g_locus_aware import LocusAwareCS2GEvaluator

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)

def main():
    logger.info("Starting Sensitivity Analysis (Window Sizes)")
    
    windows = [250, 500, 1000]
    results_summary = []
    
    for w in windows:
        logger.info(f"\nrunning window = ±{w}kb...")
        # Run Max aggregation
        evaluator = LocusAwareCS2GEvaluator(project_root, aggregation="max")
        # Ensure run() accepts window_kb
        results, summary = evaluator.run(window_kb=w)
        
        results_summary.append({
            'window_kb': w,
            'top1_acc': summary['top1_correct_mean'],
            'top1_ci_lower': summary['top1_correct_ci_lower'],
            'top1_ci_upper': summary['top1_correct_ci_upper']
        })
        
    logger.info("\nSensitivity Analysis Results:")
    df = pd.DataFrame(results_summary)
    print(df)
    
    df.to_csv(project_root / "results/cs2g_locus_aware/sensitivity_results.tsv", sep='\t', index=False)

if __name__ == "__main__":
    main()
