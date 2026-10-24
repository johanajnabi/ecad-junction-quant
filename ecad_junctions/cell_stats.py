"""Cell-level statistics with field-aware (clustered) inference. usage: python cell_stats.py results"""
import sys, numpy as np, pandas as pd
from scipy import stats
OUT = sys.argv[1] if len(sys.argv) > 1 else 'results'
c = pd.read_csv(f'{OUT}/per_cell_metrics.csv'); c['has_junction'] = c.has_junction.astype(float)
METRICS = {'junction_density': 'Junction density per cell territory (um / 100 um2)', 'has_junction': 'Fraction of cells with a detected junction', 'ridge_snr': 'Ridge SNR per cell territory',
           'continuity': 'Junction continuity (cells with junction only)', 'JC_ratio': 'Junction/cytoplasm ratio (cells with junction only)'}
rng = np.random.default_rng(11); NP, NB = 5000, 2000
def groups(tr, tm, m):
    d = c[(c.treatment == tr) & (c.time_min == tm)].dropna(subset=[m]); return [g[m].values for _, g in d.groupby('image') if len(g)]
def mean_pool(L): return np.concatenate(L).mean()
def clustered(A, B):
    """A, B: lists of per-field cell arrays. returns diff of pooled means, permutation p (swap whole fields), bootstrap CI"""
    obs = mean_pool(A) - mean_pool(B); allf = A + B; na = len(A); cnt = 0
    for _ in range(NP):
        idx = rng.permutation(len(allf)); d = mean_pool([allf[i] for i in idx[:na]]) - mean_pool([allf[i] for i in idx[na:]]); cnt += abs(d) >= abs(obs) - 1e-12
    bs = [mean_pool([A[i] for i in rng.integers(0, len(A), len(A))]) - mean_pool([B[i] for i in rng.integers(0, len(B), len(B))]) for _ in range(NB)]
    return obs, (cnt + 1) / (NP + 1), np.percentile(bs, [2.5, 97.5])
def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); n = len(p); adj = np.empty(n); run = 0
    for r, i in enumerate(o): run = max(run, (n - r) * p[i]); adj[i] = min(run, 1)
    return adj
rows = []
for m in METRICS:
    fam = []
    for tm in [0, 2, 5, 10, 15]:
        fam.append(('CalA vs Control (same time)', 'CalA', tm, 'Control', tm))
    for tr in ['Control', 'CalA']:
        for tm in [2, 5, 10, 15]: fam.append((f'{tr}: time vs own 0 min', tr, tm, tr, 0))
    out = []
    for kind, t1, tm1, t0, tm0 in fam:
        A = groups(t1, tm1, m); B = groups(t0, tm0, m)
        a = np.concatenate(A); b = np.concatenate(B)
        if m == 'has_junction':
            # chi-square-free: naive test = Mann-Whitney on 0/1 is fine; report the same
            pass
        u = stats.mannwhitneyu(a, b, alternative='two-sided').pvalue
        diff, pperm, ci = clustered(A, B)
        out.append(dict(metric=m, comparison=kind, group_A=f'{t1} {tm1} min', group_B=f'{t0} {tm0} min', n_cells_A=len(a), n_cells_B=len(b), n_fields_A=len(A), n_fields_B=len(B),
                        mean_A=a.mean(), mean_B=b.mean(), median_A=np.median(a), median_B=np.median(b), mean_diff_A_minus_B=diff, CI95_lo_cluster_bootstrap=ci[0], CI95_hi_cluster_bootstrap=ci[1],
                        p_naive_MannWhitney_cells=u, p_field_aware_permutation=pperm))
    O = pd.DataFrame(out)
    for k, g in O.groupby('comparison'):
        O.loc[g.index, 'p_naive_holm'] = holm(g.p_naive_MannWhitney_cells); O.loc[g.index, 'p_field_aware_holm'] = holm(g.p_field_aware_permutation)
    rows.append(O); print(m, flush=True)
S = pd.concat(rows); S.to_csv(f'{OUT}/cell_level_stats.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(S[S.metric.isin(['junction_density', 'has_junction'])][['metric', 'comparison', 'group_A', 'group_B', 'n_cells_A', 'n_cells_B', 'mean_A', 'mean_B', 'p_naive_holm', 'p_field_aware_holm']].round(4).to_string())
