import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os; R=os.environ.get('RES','results'); O=os.environ.get('OUTD','figures')
c=pd.read_csv(f'{R}/per_cell_metrics.csv'); os.makedirs(O,exist_ok=True); S=pd.read_csv(f'{R}/cell_level_stats.csv')
plt.rcParams.update({'font.family':['Arial','Liberation Sans','DejaVu Sans'],'font.size':13,'axes.linewidth':1.3,'xtick.major.width':1.3,'ytick.major.width':1.3})
T=[0,2,5,10,15]; NAVY='#0b0b8c'; RED='#ee1111'; PAL={'Control':NAVY,'CalA':RED}; rng=np.random.default_rng(1)
def stars(p): return 'ns' if p>=.05 else '*' if p>=.01 else '**' if p>=.001 else '***' if p>=.0001 else '****'
def make(metric,ylabel,fname,ylim,fieldlevel=False):
    fig,ax=plt.subplots(figsize=(7,5.4)); off=.2; top=ylim[1]
    for i,tm in enumerate(T):
        for tr,o in [('Control',-off),('CalA',off)]:
            d=c[(c.treatment==tr)&(c.time_min==tm)]
            if fieldlevel:
                fm=d.groupby('image')[metric].mean(); pts=fm.values; med=np.median(pts)
                ax.scatter(i+o+rng.uniform(-.14,.14,len(pts)),pts,s=22,color=PAL[tr],lw=0,zorder=2)
            else:
                v=d[metric].dropna(); pts=v.values; med=np.median(pts)
                ax.scatter(i+o+rng.uniform(-.16,.16,len(v)),np.clip(v,0,ylim[1]),s=5,color=PAL[tr],alpha=.45,lw=0,zorder=2)
                fm=d.groupby('image')[metric].median().dropna()
                ax.scatter(i+o+rng.uniform(-.05,.05,len(fm)),fm.values,s=34,facecolor='white',edgecolor='black',lw=1,zorder=4)
            ax.hlines(med,i+o-.18,i+o+.18,color='black',lw=3.2,zorder=5)
        row=S[(S.metric==metric)&(S.comparison=='CalA vs Control (same time)')&(S.group_A==f'CalA {tm} min')]
        p=row.p_field_aware_holm.iloc[0]; y=top*1.0
        ax.plot([i-off,i+off],[y,y],color='black',lw=1.2,clip_on=False)
        ax.text(i,y+.01*top,stars(p),ha='center',va='bottom',fontsize=15 if p<.05 else 12,fontweight='bold' if p<.05 else 'normal')
    m0=c[(c.treatment=='Control')&(c.time_min==0)][metric].median() if not fieldlevel else c[(c.treatment=='Control')&(c.time_min==0)].groupby('image')[metric].mean().median()
    ax.axhline(m0,color='black',ls=(0,(3,3)),lw=1.1,zorder=1)
    ax.set_xticks(range(5)); ax.set_xticklabels(T); ax.set_xlabel('Minutes post wounding'); ax.set_ylabel(ylabel); ax.set_ylim(0,ylim[1]); ax.set_xlim(-.6,4.6)
    ax.spines[['top','right']].set_visible(False)
    ax.legend(handles=[Patch(facecolor=NAVY,edgecolor='black',label='Control'),Patch(facecolor=RED,edgecolor='black',label='Calyculin A')],loc='lower center',bbox_to_anchor=(.5,1.1),ncol=2,frameon=False,handlelength=1.3)
    fig.tight_layout(); fig.savefig(fname,dpi=250); plt.close(fig)
make('junction_density','E-cad junction density per cell\n(µm / 100 µm²)',f'{O}/cells_junction_density.png',(0,20))
make('ridge_snr','E-cad ridge SNR per cell',f'{O}/cells_ridge_SNR.png',(0,3.6))
make('continuity','E-cad junction continuity per cell\n(cells with junctions only)',f'{O}/cells_continuity.png',(0,1.15))
make('has_junction','Fraction of cells with a junction\n(each dot = one field)',f'{O}/cells_fraction_with_junction.png',(0,1.15),fieldlevel=True)
print(c.junction_density.describe(), c.ridge_snr.describe())
