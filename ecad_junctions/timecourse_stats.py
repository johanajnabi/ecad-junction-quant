"""Time x treatment analysis of E-cadherin junction metrics (single biological replicate; fields = units).
usage: python timecourse_stats.py results"""
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
OUT = sys.argv[1] if len(sys.argv) > 1 else 'results'
df = pd.read_csv(f'{OUT}/per_image_metrics.csv'); df = df[df.has_wound_gap | True]
TIMES = [0, 2, 5, 10, 15]; TREATS = ['Control', 'CalA']
PRIMARY = {'skel_density': 'Junction density (um per 100 um2)', 'continuous_frac': 'Continuity (frac. in segments >=13 um)', 'ridge_snr_p95': 'Ridge SNR (95th pct)'}
SECOND = {'segments_per_100um': 'Fragmentation (segments / 100 um)', 'mean_segment_um': 'Mean segment length (um)', 'JC_ratio': 'Junction / cytoplasm ratio', 'ridge_contrast': 'Junction / local background'}
ALL = {**PRIMARY, **SECOND}
MARGIN = 0.75          # equivalence margin on the fold-change scale (+-25%): CalA 'preserved' if the 90% CI of fold-change vs its own 0 min lies in [0.75, 1/0.75]
NB = 4000; rng = np.random.default_rng(7)
def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); n = len(p); adj = np.empty(n); run = 0
    for r, i in enumerate(o): run = max(run, (n - r) * p[i]); adj[i] = min(run, 1)
    return adj
def fold_ci(a, b, level=95):
    """fold change of medians a/b with bootstrap CI"""
    a = np.asarray(a); b = np.asarray(b)
    ba = np.median(rng.choice(a, (NB, len(a))), 1); bb = np.median(rng.choice(b, (NB, len(b))), 1)
    r = ba / np.where(bb == 0, np.nan, bb); lo, hi = np.nanpercentile(r, [(100 - level) / 2, 100 - (100 - level) / 2])
    return np.median(a) / np.median(b), lo, hi
def cliffs(a, b):
    a = np.asarray(a); b = np.asarray(b); return (np.sum(a[:, None] > b[None, :]) - np.sum(a[:, None] < b[None, :])) / (len(a) * len(b))
def g(t, tm, m): return df[(df.treatment == t) & (df.time_min == tm)][m].dropna().values
rows_desc, rows_vs0, rows_trend, rows_cmp, rows_int = [], [], [], [], []
for m in ALL:
    for t in TREATS:
        base = g(t, 0, m)
        for tm in TIMES:
            x = g(t, tm, m)
            q = np.percentile(x, [25, 75])
            rows_desc.append(dict(metric=m, treatment=t, time_min=tm, n=len(x), median=np.median(x), q25=q[0], q75=q[1]))
        # per-time vs own 0 min
        ps = []; tmp = []
        for tm in TIMES[1:]:
            x = g(t, tm, m); u = stats.mannwhitneyu(x, base, alternative='two-sided')
            fc, lo, hi = fold_ci(x, base, 95); _, lo90, hi90 = fold_ci(x, base, 90)
            tmp.append(dict(metric=m, treatment=t, time_min=tm, n=len(x), n0=len(base), fold_vs_own_0min=fc, fold_CI95_lo=lo, fold_CI95_hi=hi,
                            fold_CI90_lo=lo90, fold_CI90_hi=hi90, cliffs_delta=cliffs(x, base), p_MW=u.pvalue,
                            equivalent_within_margin=bool(lo90 >= MARGIN and hi90 <= 1 / MARGIN)))
            ps.append(u.pvalue)
        for r, a in zip(tmp, holm(ps)): r['p_MW_holm_over_times'] = a
        rows_vs0 += tmp
        # monotonic trend over time (Spearman on fields)
        d = df[df.treatment == t].dropna(subset=[m]); rho = stats.spearmanr(d.time_min, d[m])
        rows_trend.append(dict(metric=m, treatment=t, n=len(d), spearman_rho_vs_time=rho.statistic, p=rho.pvalue))
    # CalA vs control at each time
    ps = []; tmp = []
    for tm in TIMES:
        a = g('CalA', tm, m); b = g('Control', tm, m); u = stats.mannwhitneyu(a, b, alternative='two-sided'); fc, lo, hi = fold_ci(a, b)
        tmp.append(dict(metric=m, time_min=tm, n_CalA=len(a), n_ctrl=len(b), median_CalA=np.median(a), median_ctrl=np.median(b), CalA_over_ctrl=fc, CI95_lo=lo, CI95_hi=hi, cliffs_delta=cliffs(a, b), p_MW=u.pvalue)); ps.append(u.pvalue)
    for r, a in zip(tmp, holm(ps)): r['p_MW_holm_over_times'] = a
    rows_cmp += tmp
    # interaction: mean over 5,10,15 min of [fold vs own 0 min], control minus CalA; permutation of treatment labels within time point
    def stat(lbl_by_time):
        out = []
        for tm in TIMES:
            pass
        return None
    def fold_stat(tbl):  # tbl: dict (treat,time)->array
        f = {}
        for t in TREATS:
            b = np.median(tbl[(t, 0)]); f[t] = np.mean([np.median(tbl[(t, tm)]) / b for tm in (5, 10, 15)])
        return f['CalA'] - f['Control']      # >0 means CalA retained more than control
    tbl = {(t, tm): g(t, tm, m) for t in TREATS for tm in TIMES}
    obs = fold_stat(tbl); cnt = 0; N = 5000
    for _ in range(N):
        tb2 = {}
        for tm in TIMES:
            pool = np.r_[tbl[('Control', tm)], tbl[('CalA', tm)]]; rng.shuffle(pool); nc = len(tbl[('Control', tm)])
            tb2[('Control', tm)] = pool[:nc]; tb2[('CalA', tm)] = pool[nc:]
        if abs(fold_stat(tb2)) >= abs(obs) - 1e-12: cnt += 1
    # bootstrap CI of the statistic
    bs = []
    for _ in range(1500):
        bs.append(fold_stat({k: rng.choice(v, len(v)) for k, v in tbl.items()}))
    rows_int.append(dict(metric=m, retained_fold_ctrl_mean_5_10_15=np.mean([np.median(tbl[('Control', tm)]) / np.median(tbl[('Control', 0)]) for tm in (5, 10, 15)]),
                         retained_fold_CalA_mean_5_10_15=np.mean([np.median(tbl[('CalA', tm)]) / np.median(tbl[('CalA', 0)]) for tm in (5, 10, 15)]),
                         difference_CalA_minus_ctrl=obs, CI95_lo=np.percentile(bs, 2.5), CI95_hi=np.percentile(bs, 97.5), p_perm_interaction=(cnt + 1) / (N + 1)))
D = pd.DataFrame(rows_desc); V = pd.DataFrame(rows_vs0); T = pd.DataFrame(rows_trend); C = pd.DataFrame(rows_cmp); I = pd.DataFrame(rows_int)
I['p_interaction_holm_primary'] = np.nan; pm = I.metric.isin(PRIMARY)
I.loc[pm, 'p_interaction_holm_primary'] = holm(I.loc[pm, 'p_perm_interaction'])
# voltage-adjusted sensitivity: residualise log(metric) on E-cad PMT voltage within treatment, then trend vs time
sens = []
for m in ALL:
    for t in TREATS:
        d = df[df.treatment == t].dropna(subset=[m, 'ecad_pmt_V']); y = np.log(np.clip(d[m], 1e-6, None))
        X = np.c_[np.ones(len(d)), d.ecad_pmt_V]; b = np.linalg.lstsq(X, y, rcond=None)[0]; res = y - X @ b
        r1 = stats.spearmanr(d.time_min, d[m]); r2 = stats.spearmanr(d.time_min, res)
        sens.append(dict(metric=m, treatment=t, rho_raw=r1.statistic, p_raw=r1.pvalue, rho_after_PMT_adjust=r2.statistic, p_after_PMT_adjust=r2.pvalue,
                         rho_PMT_vs_time=stats.spearmanr(d.time_min, d.ecad_pmt_V).statistic))
S = pd.DataFrame(sens)
with pd.ExcelWriter(f'{OUT}/timecourse_stats.xlsx') as xw:
    for n, t in [('descriptives', D), ('vs_own_0min', V), ('trend_over_time', T), ('CalA_vs_control', C), ('interaction', I), ('PMT_sensitivity', S)]: t.round(5).to_excel(xw, sheet_name=n, index=False)
for n, t in [('descriptives', D), ('vs_own_0min', V), ('trend_over_time', T), ('CalA_vs_control', C), ('interaction', I), ('PMT_sensitivity', S)]: t.to_csv(f'{OUT}/stats_{n}.csv', index=False)
# ---------- figures ----------
COL = {'Control': '#2a78d6', 'CalA': '#eb6834'}; OFF = {'Control': -0.18, 'CalA': 0.18}
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#52514e', 'axes.labelcolor': '#0b0b0b', 'xtick.color': '#52514e', 'ytick.color': '#52514e'})
xs = {tm: i for i, tm in enumerate(TIMES)}
def tc_panel(ax, m, lab, fold=False):
    for t in TREATS:
        meds = []; los = []; his = []
        for tm in TIMES:
            x = g(t, tm, m); base = np.median(g(t, 0, m)) if fold else 1.0; x = x / base
            jit = rng.normal(0, 0.035, len(x)); ax.scatter(xs[tm] + OFF[t] + jit, x, s=9, color=COL[t], alpha=.35, linewidths=0, zorder=2)
            meds.append(np.median(x)); los.append(np.percentile(x, 25)); his.append(np.percentile(x, 75))
        px = [xs[tm] + OFF[t] for tm in TIMES]
        ax.plot(px, meds, '-o', color=COL[t], lw=2, ms=6, markeredgecolor='white', markeredgewidth=1.2, zorder=4, label='Calyculin A' if t == 'CalA' else 'Control')
        ax.vlines(px, los, his, color=COL[t], lw=2, zorder=3)
    if fold: ax.axhline(1, color='#9a9992', lw=.8, ls=':'); ax.axhspan(MARGIN, 1 / MARGIN, color='#9a9992', alpha=.12, lw=0)
    ax.set_xticks(range(len(TIMES))); ax.set_xticklabels([str(t) for t in TIMES]); ax.set_xlabel('Minutes after wounding'); ax.set_title(lab, fontsize=9.5, loc='left')
    if fold: ax.set_ylabel('Fold vs own 0 min')
for kind, fold in [('abs', False), ('fold', True)]:
    fig, axs = plt.subplots(2, 4, figsize=(15.5, 7.2)); axs = axs.ravel()
    for ax, (m, l) in zip(axs, ALL.items()): tc_panel(ax, m, l, fold)
    h, lb = axs[0].get_legend_handles_labels(); axs[-1].axis('off'); axs[-1].legend(h, lb, loc='center', frameon=False, fontsize=11, title='Median and IQR; dots = fields')
    fig.suptitle('E-cadherin junction metrics across the wound time course' + (' (each treatment normalised to its own 0 min median; grey band = +-25%)' if fold else ''), fontsize=11, x=0.01, ha='left')
    fig.tight_layout(); fig.savefig(f'{OUT}/fig_timecourse_{kind}.png', dpi=150); plt.close(fig)
# voltage confound figure
fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
for t in TREATS:
    d = df[df.treatment == t]
    for tm in TIMES:
        x = d[d.time_min == tm].ecad_pmt_V; axs[0].scatter(np.full(len(x), xs[tm]) + OFF[t] + rng.normal(0, .03, len(x)), x, s=14, color=COL[t], alpha=.6, linewidths=0)
    axs[1].scatter(d.ecad_pmt_V, d.skel_density, s=14, color=COL[t], alpha=.6, linewidths=0, label='Calyculin A' if t == 'CalA' else 'Control')
axs[0].set_xticks(range(5)); axs[0].set_xticklabels([str(t) for t in TIMES]); axs[0].set_xlabel('Minutes after wounding'); axs[0].set_ylabel('E-cad detector voltage (V)'); axs[0].set_title('Acquisition voltage by group', loc='left', fontsize=9.5)
axs[1].set_xlabel('E-cad detector voltage (V)'); axs[1].set_ylabel('Junction density'); axs[1].legend(frameon=False); axs[1].set_title('Voltage vs junction density', loc='left', fontsize=9.5)
fig.tight_layout(); fig.savefig(f'{OUT}/fig_voltage_confound.png', dpi=150)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print('N per group:\n', df.groupby(['treatment', 'time_min']).size().unstack())
print(D[D.metric == 'skel_density'].round(3).to_string()); print(V[V.metric.isin(PRIMARY)].round(4).to_string()); print(T.round(4).to_string()); print(I.round(4).to_string()); print(S[S.metric.isin(PRIMARY)].round(3).to_string())
