import sys; sys.path.insert(0,'.')
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import style as S; from style import *
S.apply()

X = ['Receptor-\ngrouped', '50%', '40%', '30%']
series = {
 'Full sequence':        (SEQ,   [0.790,0.638,0.612,0.539], [0.744,0.584,0.560,0.480], [0.832,0.690,0.663,0.600], '-'),
 'Interface sequence':   ('#B07CD6',[0.779,0.649,0.610,0.552],[0.729,0.591,0.560,0.502],[0.827,0.703,0.664,0.606], (0,(4,1.6))),
 'Endpoint geometry':    (GEOM,  [0.776,0.731,0.683,0.517], [0.726,0.676,0.627,0.456], [0.826,0.784,0.737,0.578], '-'),
 'Global identity kNN':  ('#8D99AE',[0.802,0.570,0.552,0.550],[0.756,0.515,0.500,0.494],[0.844,0.623,0.602,0.605], (0,(1.6,1.6))),
}

fig = plt.figure(figsize=(7.0, 2.95))
gs = GridSpec(1, 3, figure=fig, width_ratios=[1.30, 0.95, 1.05], wspace=0.40,
              left=0.075, right=0.995, top=0.86, bottom=0.235)

# ---------------- A ----------------
axA = fig.add_subplot(gs[0,0])
x = np.arange(4)
for name, (c, v, l, h, ls) in series.items():
    axA.fill_between(x, l, h, color=c, alpha=0.13, lw=0, zorder=2)
    axA.plot(x, v, color=c, lw=1.5, ls=ls, marker='o', ms=3.6,
             mec='white', mew=0.8, zorder=4, label=name)
axA.axhline(0.5, color='#9AA0A6', lw=0.8, ls=(0,(4,2)), zorder=1)
axA.text(2.15, 0.508, 'chance', fontsize=5.8, color='#6C757D', ha='center')
axA.set_xticks(x); axA.set_xticklabels(X, fontsize=6.4, linespacing=1.15)
axA.set_xlim(-0.28, 3.45); axA.set_ylim(0.44, 0.90)
axA.set_ylabel('Macro ROC AUC')
axA.set_xlabel('Sequence-identity separation between folds', labelpad=2)
axA.legend(frameon=False, loc='lower left', fontsize=5.9, handlelength=1.6,
           handletextpad=0.4, labelspacing=0.22, borderpad=0.0)
axA.annotate('geometry > sequence', xy=(1.06, 0.742), xytext=(1.35, 0.845),
             fontsize=6.0, color=GEOM, fontweight='bold',
             arrowprops=dict(arrowstyle='-', lw=0.7, color=GEOM))
S.hgrid(axA); S.panel(axA, 'A', x=-0.20)

# ---------------- B ----------------
axB = fig.add_subplot(gs[0,1])
d  = [0.014, -0.093, -0.070, 0.022]
dl = [-0.037, -0.162, -0.141, -0.061]
dh = [0.065, -0.025, 0.004, 0.101]
cols = [SEQ if v > 0 else GEOM for v in d]
axB.axhline(0, color='#C1121F', lw=1.0, zorder=3)
axB.bar(x, d, width=0.58, color=cols, edgecolor='white', lw=0.7, zorder=4)
axB.errorbar(x, d, yerr=[np.array(d)-np.array(dl), np.array(dh)-np.array(d)],
             fmt='none', ecolor=INK, elinewidth=0.9, capsize=2.2, capthick=0.9, zorder=5)
axB.set_xticks(x); axB.set_xticklabels(X, fontsize=6.4, linespacing=1.15)
axB.set_xlim(-0.55, 3.55); axB.set_ylim(-0.20, 0.14)
axB.set_ylabel(r'$\Delta$ AUC (sequence $-$ geometry)', fontsize=6.9)
axB.text(0.03, 0.955, 'sequence better', transform=axB.transAxes, fontsize=5.9,
         color=SEQ, fontweight='bold', va='top')
axB.text(0.03, 0.045, 'geometry better', transform=axB.transAxes, fontsize=5.9,
         color=GEOM, fontweight='bold', va='bottom')
S.hgrid(axB); S.panel(axB, 'B', x=-0.30)

# ---------------- C ----------------
axC = fig.add_subplot(gs[0,2])
lab = ['Interface','Matched\nnoninterface','Global\nidentity','Endpoint\ngeometry']
v  = [0.552, 0.556, 0.550, 0.517]
lo = [0.502, 0.500, 0.494, 0.456]
hi = [0.606, 0.613, 0.605, 0.578]
cc = ['#B07CD6', '#C9C2D6', '#8D99AE', GEOM]
xx = np.arange(4)
axC.bar(xx, np.array(v)-0.44, bottom=0.44, width=0.60, color=cc,
        edgecolor='white', lw=0.7, zorder=3)
axC.errorbar(xx, v, yerr=[np.array(v)-np.array(lo), np.array(hi)-np.array(v)],
             fmt='none', ecolor=INK, elinewidth=0.9, capsize=2.2, capthick=0.9, zorder=4)
axC.axhline(0.5, color='#9AA0A6', lw=0.9, ls=(0,(4,2)), zorder=2)
axC.set_xticks(xx); axC.set_xticklabels(lab, fontsize=5.9, linespacing=1.1)
axC.set_ylim(0.44, 0.68)
axC.set_yticks([0.45,0.50,0.55,0.60,0.65]); axC.set_xlim(-0.62, 3.62)
axC.set_ylabel('Macro ROC AUC', fontsize=6.9)
axC.set_title('At 30% identity separation', fontsize=6.8, color=INK, pad=13)
axC.text(1.5, 0.663, r'interface $-$ noninterface:', ha='center', fontsize=5.9, color='#495057')
axC.text(1.5, 0.645, r'$-$0.004 ($-$0.068 to 0.062)', ha='center', fontsize=5.9,
         color='#495057', fontweight='bold')
S.hgrid(axC); S.panel(axC, 'C', x=-0.28)

for ext in ('pdf','svg','png'):
    fig.savefig(f'../figures/Figure_3.{ext}', bbox_inches='tight', facecolor='white')
print('fig3 done')
