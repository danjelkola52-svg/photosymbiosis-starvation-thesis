"""Thesis physiology figures (Figures 4-11; Figure 3 is the phylogeny): one consistent style, no panel titles,
Holm-adjusted Welch p-values printed to three decimals.

Usage:  python thesis_physiology_figures.py <master_data_all_parameters_1.xlsx> <output folder>
"""
import sys, os, warnings
import numpy as np, pandas as pd, matplotlib
from scipy import stats
warnings.filterwarnings('ignore')
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

SRC = sys.argv[1] if len(sys.argv) > 1 else 'master_data_all_parameters_1.xlsx'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'thesis_figures'
os.makedirs(OUT, exist_ok=True)

for fam in ['Arial', 'Liberation Sans', 'DejaVu Sans']:
    if any(fam == f.name for f in matplotlib.font_manager.fontManager.ttflist):
        plt.rcParams['font.family'] = fam; break
plt.rcParams.update({'mathtext.default': 'regular', 'axes.linewidth': 0.8})

PAL = {'Fed': '#5B8FC7', 'Starved': '#E08A62'}
HOSTS = ['Exaiptasia', 'Waminoa']

d = pd.read_excel(SRC, sheet_name='Master')
d['Organism'] = d['Organism'].replace({'Aiptasia': 'Exaiptasia'})

# (figure number, file stub, column, y-axis label, scale factor)
FIGS = [
    (4,  'Fig04_TotalProtein',       'Protein (μg / animal)',            'Total protein (µg animal$^{-1}$)', 1),
    (5,  'Fig05_SymbiontsPerAnimal', 'Symbiont cells (per animal)',      'Symbionts per animal (× 10$^{3}$)', 1e-3),
    (6,  'Fig06_SymbiontDensity',    'Cells / μg protein',               'Symbiont density\n(cells µg$^{-1}$ protein)', 1),
    (7,  'Fig07_NetO2',              'ΔO₂ absolute (µg O₂ h⁻¹)',         'Net O$_2$ production\n(µg O$_2$ h$^{-1}$ animal$^{-1}$)', 1),
    (8,  'Fig08_O2perProtein',       'ΔO₂ (nmol O₂ h⁻¹ µg protein⁻¹)',   'Net O$_2$ per protein\n(nmol O$_2$ h$^{-1}$ µg$^{-1}$)', 1),
    (11, 'Fig11_O2perCell',          'ΔO₂ (nmol O₂ h⁻¹ cell⁻¹)',         'Net O$_2$ per symbiont cell\n(pmol O$_2$ h$^{-1}$ cell$^{-1}$)', 1e3),
    (9,  'Fig09_FvFm',               'PAM Fv/Fm',                        'Maximum quantum yield (F$_v$/F$_m$)', 1),
    (10, 'Fig10_NH4perProtein',      'NH₄⁺ / μg protein (μM μg⁻¹)',      'NH$_4^+$ uptake (µM µg$^{-1}$ protein)', 1),
]

def fmt_p(p):
    return 'p < 0.001' if p < 0.001 else f'p = {p:.3f}'

def holm(ps):
    order = np.argsort(ps); adj = np.empty(len(ps)); run = 0
    for k, i in enumerate(order):
        run = max(run, min(1, ps[i] * (len(ps) - k))); adj[i] = run
    return adj

rows = []
rng = np.random.default_rng(7)
for num, stub, col, ylab, scale in FIGS:
    data = {h: {t: d[(d.Organism == h) & (d.Treatment == t)][col].dropna().values * scale for t in ['Fed', 'Starved']} for h in HOSTS}
    praw = [stats.ttest_ind(data[h]['Fed'], data[h]['Starved'], equal_var=False).pvalue for h in HOSTS]
    padj = holm(np.array(praw))
    allv = np.concatenate([np.concatenate(list(data[h].values())) for h in HOSTS])
    lo, hi = allv.min(), allv.max(); span = hi - lo
    ylim = (lo - 0.08 * span, hi + 0.26 * span)

    fig, axes = plt.subplots(1, 2, figsize=(6.3, 3.4), sharey=True)
    for k, (ax, h) in enumerate(zip(axes, HOSTS)):
        for x, t in enumerate(['Fed', 'Starved']):
            v = data[h][t]
            ax.boxplot(v, positions=[x], widths=0.55, patch_artist=True, showfliers=False,
                       boxprops=dict(facecolor=PAL[t], alpha=0.85, edgecolor='black', linewidth=0.7),
                       medianprops=dict(color='black', linewidth=1.3),
                       whiskerprops=dict(color='black', linewidth=0.7), capprops=dict(color='black', linewidth=0.7))
            jx = x + rng.uniform(-0.10, 0.10, len(v))
            ax.scatter(jx, v, s=22, c=PAL[t], edgecolors='black', linewidths=0.5, zorder=3)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f'Fed\n(n = {len(data[h]["Fed"])})', f'Starved\n(n = {len(data[h]["Starved"])})'], fontsize=8.5)
        ax.set_xlim(-0.6, 1.6); ax.set_ylim(*ylim)
        yb = ylim[1] - 0.10 * (ylim[1] - ylim[0]); tick = 0.022 * (ylim[1] - ylim[0])
        ax.plot([0, 0, 1, 1], [yb - tick, yb, yb, yb - tick], color='black', lw=0.7)
        p = padj[k]
        ax.text(0.5, yb + 0.02 * (ylim[1] - ylim[0]), fmt_p(p), ha='center', va='bottom', fontsize=8.5,
                fontweight='bold' if p < 0.05 else 'normal')
        for sp in ['top', 'right']: ax.spines[sp].set_visible(False)
        ax.tick_params(labelsize=8.5)
        ax.yaxis.grid(True, color='0.92', lw=0.5); ax.set_axisbelow(True)
        ax.text(-0.16 if k == 0 else -0.06, 1.02, 'AB'[k], transform=ax.transAxes, fontsize=12, fontweight='bold')
        rows.append(dict(Figure=num, Parameter=stub.split('_', 1)[1], Host=h, n_Fed=len(data[h]['Fed']),
                         n_Starved=len(data[h]['Starved']), p_Welch=praw[k], p_Holm=p))
    axes[0].set_ylabel(ylab, fontsize=9)
    if max(abs(np.array(ylim))) >= 1e4:
        axes[0].yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:,.0f}'))
    fig.tight_layout(w_pad=2.0)
    fig.savefig(os.path.join(OUT, stub + '.png'), dpi=600)
    fig.savefig(os.path.join(OUT, stub + '.pdf'))
    plt.close(fig)

pd.DataFrame(rows).to_csv(os.path.join(OUT, 'figure_pvalues.csv'), index=False)
print('done:', OUT)
