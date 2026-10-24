"""Prism/Fig-2E-style swarm plot: Control vs Calyculin A across minutes after wounding. Usage: python plot_fig2E_style.py results"""
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import seaborn as sns
OUT = sys.argv[1] if len(sys.argv) > 1 else 'results'
df = pd.read_csv(f'{OUT}/per_image_metrics.csv'); C = pd.read_csv(f'{OUT}/stats_CalA_vs_control.csv')
plt.rcParams.update({'font.family': ['Arial', 'Liberation Sans', 'DejaVu Sans'], 'font.size': 13, 'axes.linewidth': 1.3, 'xtick.major.width': 1.3, 'ytick.major.width': 1.3,
                     'xtick.direction': 'out', 'ytick.direction': 'out'})
TIMES = [0, 2, 5, 10, 15]; NAVY = '#0b0b8c'; RED = '#ee1111'
PAL = {'Control': NAVY, 'CalA': RED}; LAB = {'Control': 'Control', 'CalA': 'Calyculin A'}
def stars(p): return 'ns' if p >= 0.05 else '*' if p >= 0.01 else '**' if p >= 0.001 else '***' if p >= 0.0001 else '****'
def make(metric, ylabel, fname, ylim=None, ref=True):
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    d = df.copy(); d['time_min'] = pd.Categorical(d.time_min, TIMES, ordered=True)
    sns.swarmplot(data=d, x='time_min', y=metric, hue='treatment', hue_order=['Control', 'CalA'], palette=PAL, dodge=True, size=3.3, ax=ax, linewidth=0, zorder=2)
    ax.legend_.remove()
    off = 0.2; top = d[metric].max() if ylim is None else ylim[1]
    for i, tm in enumerate(TIMES):
        for t, o in [('Control', -off), ('CalA', off)]:
            x = d[(d.treatment == t) & (d.time_min == tm)][metric].dropna()
            ax.hlines(x.median(), i + o - 0.17, i + o + 0.17, color='black', lw=3.2, zorder=5)
        p = C[(C.metric == metric) & (C.time_min == tm)].p_MW_holm_over_times.iloc[0]
        y = d[d.time_min == tm][metric].max() * 1.05 + 0.02 * top
        ax.plot([i - off, i + off], [y, y], color='black', lw=1.2, clip_on=False)
        ax.text(i, y + 0.012 * top, stars(p), ha='center', va='bottom', fontsize=15 if p < 0.05 else 12, fontweight='bold' if p < 0.05 else 'normal')
    if ref:
        m0 = d[(d.treatment == 'Control') & (d.time_min == 0)][metric].median(); ax.axhline(m0, color='black', ls=(0, (3, 3)), lw=1.1, zorder=1)
    ax.set_xlabel('Minutes post wounding', labelpad=6); ax.set_ylabel(ylabel)
    ax.set_xticks(range(len(TIMES))); ax.set_xticklabels([str(t) for t in TIMES])
    ax.set_ylim(0, ylim[1] if ylim else None); ax.spines[['top', 'right']].set_visible(False)
    # legend (coloured squares) above the plot
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=NAVY, edgecolor='black', label='Control'), Patch(facecolor=RED, edgecolor='black', label='Calyculin A')], loc='lower center', bbox_to_anchor=(0.5, 1.12), ncol=2, frameon=False, handlelength=1.3, columnspacing=1.6)
    fig.tight_layout(); fig.savefig(fname, dpi=300); plt.close(fig)
make('skel_density', 'E-cad junction density (µm / 100 µm²)', f'{OUT}/fig2E_style_junction_density.png', ylim=(0, 18))
make('continuous_frac', 'E-cad junction continuity (fraction)', f'{OUT}/fig2E_style_continuity.png', ylim=(0, 1.1))
make('ridge_snr_p95', 'E-cad ridge SNR (95th pct)', f'{OUT}/fig2E_style_ridge_SNR.png', ylim=(0, 6.5))
