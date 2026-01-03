
import pandas as pd
import numpy as np
from scipy.stats import binomtest

def run_mcnemar():
    df = pd.read_csv('results/baselines/post2021_predictions_all_methods.tsv', sep='\t')
    
    methods = df['method'].unique()
    print(f"Methods: {methods}")
    
    # Pivot to get Correct/Incorrect matrix
    wide = df.pivot(index='locus_id', columns='method', values='top1_correct')
    
    target = 'cS2G_LocusAware_max'
    if target not in wide.columns:
        print(f"Method {target} not found in results.")
        return

    print(f"\nComparing {target} vs others (McNemar's test):")
    print(f"{'Comparator':<20} | {'Both 1':<6} | {'Only A':<6} | {'Only B':<6} | {'Both 0':<6} | {'p-value':<10} | {'Significance'}")
    print("-" * 90)
    
    for other in wide.columns:
        if other == target:
            continue
            
        # cell [0,1]: Target wrong, Other correct (Only B) => 'c'
        # cell [1,0]: Target correct, Other wrong (Only A) => 'b'
        
        only_other = ((wide[target] == 0) & (wide[other] == 1)).sum() # Disagreements where Other is right
        only_target = ((wide[target] == 1) & (wide[other] == 0)).sum() # Disagreements where Target is right
        both_correct = ((wide[target] == 1) & (wide[other] == 1)).sum()
        both_wrong = ((wide[target] == 0) & (wide[other] == 0)).sum()
        
        # McNemar Exact Test (Binomial test on discordant pairs)
        # H0: p(only_target) = p(only_other) = 0.5
        n_discordant = only_target + only_other
        
        if n_discordant == 0:
            p_val = 1.0
        else:
            # binomtest(k, n, p=0.5)
            # We want two-sided
            res = binomtest(only_target, n_discordant, 0.5, alternative='two-sided')
            p_val = res.pvalue
            
        sig = "**" if p_val < 0.005 else "*" if p_val < 0.05 else "ns"
        
        print(f"{other:<20} | {both_correct:<6} | {only_target:<6} | {only_other:<6} | {both_wrong:<6} | {p_val:.2e}   | {sig}")

if __name__ == "__main__":
    run_mcnemar()
