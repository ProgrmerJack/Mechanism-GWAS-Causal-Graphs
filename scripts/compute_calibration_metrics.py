"""Compute L2G v25.12 vs v22.09 calibration metrics"""
import pandas as pd
import numpy as np
from sklearn.calibration import calibration_curve
from sklearn.metrics import brier_score_loss
from sklearn.linear_model import LogisticRegression

# Load matched data
df_merged = pd.read_csv('validation_bundle/calibration/gene_predictions_with_l2g_versions.tsv', sep='\t')
df_complete = df_merged.dropna(subset=['calibrated_prob', 'l2g_v22_max', 'l2g_v25_max', 'truth'])
print(f'Complete cases: {len(df_complete):,}')

# Fixed ECE computation using quantile strategy
def compute_ece(y_true, y_prob, n_bins=10):
    prob_true, prob_pred = calibration_curve(y_true, y_prob, n_bins=n_bins, strategy='quantile')
    ece = np.mean(np.abs(prob_true - prob_pred))
    return ece

results = {}

# Mechanism Graphs
print("\nComputing Mechanism Graphs metrics...")
results['Mechanism Graphs'] = {
    'ECE': compute_ece(df_complete['truth'], df_complete['calibrated_prob']),
    'Brier': brier_score_loss(df_complete['truth'], df_complete['calibrated_prob']),
    'n': len(df_complete)
}
X = df_complete['calibrated_prob'].values.reshape(-1, 1)
model = LogisticRegression(penalty=None, max_iter=1000)
model.fit(X, df_complete['truth'])
results['Mechanism Graphs']['Slope'] = model.coef_[0][0]

# L2G v22.09
print("Computing L2G v22.09 metrics...")
results['L2G v22.09'] = {
    'ECE': compute_ece(df_complete['truth'], df_complete['l2g_v22_max']),
    'Brier': brier_score_loss(df_complete['truth'], df_complete['l2g_v22_max']),
    'n': len(df_complete)
}
X = df_complete['l2g_v22_max'].values.reshape(-1, 1)
model = LogisticRegression(penalty=None, max_iter=1000)
model.fit(X, df_complete['truth'])
results['L2G v22.09']['Slope'] = model.coef_[0][0]

# L2G v25.12
print("Computing L2G v25.12 metrics...")
results['L2G v25.12'] = {
    'ECE': compute_ece(df_complete['truth'], df_complete['l2g_v25_max']),
    'Brier': brier_score_loss(df_complete['truth'], df_complete['l2g_v25_max']),
    'n': len(df_complete)
}
X = df_complete['l2g_v25_max'].values.reshape(-1, 1)
model = LogisticRegression(penalty=None, max_iter=1000)
model.fit(X, df_complete['truth'])
results['L2G v25.12']['Slope'] = model.coef_[0][0]

print('\n' + '='*80)
print('CALIBRATION METRICS COMPARISON')
print('='*80)
print(f"{'Method':<20} {'ECE':>10} {'Brier':>10} {'Slope':>10} {'n':>10}")
print('-'*80)
for method, metrics in results.items():
    print(f"{method:<20} {metrics['ECE']:>10.4f} {metrics['Brier']:>10.4f} {metrics['Slope']:>10.2f} {metrics['n']:>10,}")

print('\n' + '='*80)
print('VERSION IMPROVEMENT ANALYSIS')
print('='*80)
ece_improvement = results['L2G v22.09']['ECE'] - results['L2G v25.12']['ECE']
ece_improvement_pct = (ece_improvement / results['L2G v22.09']['ECE']) * 100
print(f"L2G ECE improvement (v22→v25): {ece_improvement:.4f} ({ece_improvement_pct:+.1f}%)")

mechanism_advantage = results['L2G v25.12']['ECE'] / results['Mechanism Graphs']['ECE']
print(f"Mechanism Graphs advantage over L2G v25: {mechanism_advantage:.2f}x better ECE")

# Save results
results_df = pd.DataFrame(results).T
results_df.to_csv('validation_bundle/calibration/l2g_version_comparison.tsv', sep='\t')
print('\nSaved results to validation_bundle/calibration/l2g_version_comparison.tsv')
