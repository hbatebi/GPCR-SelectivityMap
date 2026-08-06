import sys; sys.path.insert(0,'.')
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.gridspec import GridSpec
import style as S; from style import *
import receptor as R
S.apply()

loro = pd.read_csv('../source_data/Figure_1E_leave_one_receptor_out.csv')
mu   = pd.read_csv('../source_data/Figure_1D_activation_matching.csv')

fig = plt.figure(figsize=(7.0, 5.95))
gs = GridSpec(2, 10, figure=fig, height_ratios=[1.24, 1.0],
              hspace=0.34, wspace=2.6, left=0.075, right=0.985, top=0.93, bottom=0.085)

# ---------------- A ----------------
axA = fig.add_subplot(gs[0,0:4]); axA.set_axis_off()
import matplotlib.image as mpimg
img = mpimg.imread('../figures/Figure_1A_schematic_composite.png')
axA.imshow(img, interpolation='lanczos', aspect='equal', zorder=1)
h, w = img.shape[0], img.shape[1]
axA.set_xlim(-w*0.02, w*1.02); axA.set_ylim(h*1.15, -h*0.24)

LC, RC = 0.254*w, 0.771*w
axA.text(LC, -h*0.155, 'G protein-bound', ha='center', va='bottom',
         fontsize=7.2, fontweight='bold', color=BOUND)
axA.text(RC, -h*0.155, 'Receptor only', ha='center', va='bottom',
         fontsize=7.2, fontweight='bold', color=FREE)
axA.text(LC, -h*0.025, '15 units, 12 receptors', ha='center', va='bottom',
         fontsize=6.0, color='#495057')
axA.text(RC, -h*0.090, 'strict entity audit,', ha='center', va='bottom',
         fontsize=6.0, color='#495057')
axA.text(RC, -h*0.025, 'label-blind reference', ha='center', va='bottom',
         fontsize=6.0, color='#495057')
axA.add_patch(FancyArrowPatch((LC, h*1.035), (RC, h*1.035), arrowstyle='<|-|>',
              mutation_scale=7, lw=0.9, color='#6C757D', zorder=8))
axA.text((LC+RC)/2, h*1.065, 'matched by receptor', ha='center', va='top',
         fontsize=6.3, color='#495057')
S.panel(axA, 'A', x=-0.04, y=1.00)

# ---------------- B ----------------
axB = fig.add_subplot(gs[0,4:7])
labels = ['Bound', 'Receptor' + '\n' + 'only', 'Deform.' + '\n' + 'vector',
          'Free +' + '\n' + 'deform.']
vals = np.array([0.7952, 0.6187, 0.5680, 0.6367])
lo   = np.array([0.6546, 0.4144, 0.4131, 0.4365])
hi   = np.array([0.9444, 0.7915, 0.7062, 0.8784])
cols = [BOUND, FREE, DEFORM, COMB]
x = np.arange(4)
axB.bar(x, vals, width=0.62, color=cols, edgecolor='white', lw=0.8, zorder=3)
axB.errorbar(x, vals, yerr=[vals-lo, hi-vals], fmt='none', ecolor=INK,
             elinewidth=0.9, capsize=2.4, capthick=0.9, zorder=4)
axB.axhline(0.5, color='#9AA0A6', lw=0.8, ls=(0,(4,2)), zorder=2)
axB.text(-0.47, 0.512, 'chance', fontsize=5.5, color='#6C757D', ha='left')
for xi, v, h in zip(x, vals, hi):
    axB.text(xi, h+0.012, f'{v:.3f}', ha='center', va='bottom', fontsize=6.3,
             fontweight='bold', color=INK, zorder=5)
axB.set_xticks(x); axB.set_xticklabels(labels, fontsize=6.3, linespacing=1.15)
axB.set_ylim(0.40, 1.10); axB.set_yticks([0.4,0.5,0.6,0.7,0.8,0.9,1.0])
axB.set_ylabel('Macro ROC AUC')
S.hgrid(axB)
axB.plot([0,1],[1.045,1.045], color=INK, lw=0.8, clip_on=False)
axB.text(0.5, 1.058, r'$\Delta$ = 0.176 (0.025-0.412)', ha='center', fontsize=6.3,
         fontweight='bold', color=INK)
S.panel(axB, 'B')

# ---------------- C : per-class ----------------
axC = fig.add_subplot(gs[0,7:10])
pc = [('Gs', 5, 0.610, 0.920, GS), ('Gi/o', 6, 0.667, 0.852, GI), ('Gq/11', 4, 0.580, 0.614, GQ)]
for nm, n, f, bnd, c in pc:
    axC.plot([0, 1], [f, bnd], color=c, lw=1.6, zorder=3)
    axC.scatter([0, 1], [f, bnd], s=32, c=[c], edgecolor='white', lw=0.9, zorder=5)
    axC.text(1.07, bnd, f'{nm} (n={n})', va='center', ha='left', fontsize=6.0,
             color=c, fontweight='bold')
axC.axhline(0.5, color='#9AA0A6', lw=0.8, ls=(0,(4,2)), zorder=2)
axC.set_xticks([0, 1]); axC.set_xticklabels(['Receptor\nonly', 'Bound'], fontsize=6.2,
                                            linespacing=1.15)
axC.set_xlim(-0.28, 1.95); axC.set_ylim(0.46, 1.0)
axC.set_ylabel('One-vs-rest ROC AUC', fontsize=6.9)
axC.text(0.5, 0.965, 'gain is class-uneven', ha='center', fontsize=6.0,
         color='#495057', style='italic')
S.hgrid(axC); S.panel(axC, 'C', x=-0.30)

# ---------------- D : activation matching ----------------
axD = fig.add_subplot(gs[1,0:5])
cmap = {'Gs':GS, 'Gi/o':GI, 'Gq/11':GQ}
axD.plot([0.05,0.72],[0.05,0.72], color='#9AA0A6', lw=0.9, ls=(0,(4,2)), zorder=2)
for cls, g in mu.groupby('transducer_family'):
    axD.scatter(g['bound__activation_consensus_score'],
                g['receptor_only__activation_consensus_score'],
                s=34, c=cmap[cls], edgecolor='white', lw=0.8, zorder=4, label=cls)
for _, r_ in mu.iterrows():
    if abs(r_['bound__activation_consensus_score'] -
           r_['receptor_only__activation_consensus_score']) > 0.15:
        axD.annotate(r_['receptor_name'].split('_')[0].upper(),
                     (r_['bound__activation_consensus_score'],
                      r_['receptor_only__activation_consensus_score']),
                     textcoords='offset points', xytext=(-6,-9),
                     fontsize=5.8, color='#495057', ha='center')
axD.set_xlim(0.05,0.72); axD.set_ylim(0.05,0.72)
axD.set_xlabel('Bound activation score')
axD.set_ylabel('Receptor-only activation score')
axD.legend(frameon=False, loc='upper left', handletextpad=0.3, borderpad=0.1,
           labelspacing=0.25, fontsize=6.2)
axD.text(0.975, 0.045,
         'references are active-state entries;' + '\n' +
         'labelled points are the two least' + '\n' + 'closely matched pairs',
         transform=axD.transAxes, ha='right', va='bottom', fontsize=5.9,
         color='#495057', linespacing=1.3)
S.vgrid(axD); S.hgrid(axD)
S.panel(axD, 'D')

# ---------------- E : LORO ----------------
axE = fig.add_subplot(gs[1,5:10])
lo_ = loro.sort_values('delta')
y = np.arange(len(lo_))
axE.axvline(0, color='#C1121F', lw=1.0, zorder=3)
axE.axvline(0.176, color=INK, lw=0.9, ls=(0,(4,2)), zorder=3)
axE.scatter(lo_['delta'], y, s=30, c=BOUND, edgecolor='white', lw=0.8, zorder=5)
axE.set_yticks(y)
axE.set_yticklabels([n.split('_')[0].upper() for n in lo_['receptor']], fontsize=6.0)
axE.set_xlim(-0.025, 0.30); axE.set_ylim(-1.5, len(lo_)-0.1)
axE.set_xticks([0,0.05,0.10,0.15,0.20,0.25,0.30])
axE.set_xlabel(r'$\Delta$ macro ROC AUC (bound $-$ receptor-only)')
axE.text(0.176, len(lo_)-0.28, 'all-receptor' + '\n' + 'estimate', ha='center',
         va='top', fontsize=5.9, color=INK, linespacing=1.15)
axE.text(0.005, -0.75, 'no effect', ha='left', va='center', fontsize=5.9, color='#C1121F')
axE.text(0.985, 0.055, 'omitting each receptor in turn:' + '\n' + 'range 0.130-0.250',
         transform=axE.transAxes, ha='right', va='bottom', fontsize=6.1,
         color='#495057', fontweight='bold', linespacing=1.25)
S.vgrid(axE)
S.panel(axE, 'E')

for ext in ('pdf','svg','png'):
    fig.savefig(f'../figures/Figure_1.{ext}', bbox_inches='tight', facecolor='white')
print('fig1 done')
