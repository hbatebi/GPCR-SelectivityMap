import sys; sys.path.insert(0,'.')
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import style as S; from style import *
S.apply()

fig = plt.figure(figsize=(7.0, 3.30))
gs = GridSpec(1, 3, figure=fig, width_ratios=[1.10, 1.16, 1.02], wspace=0.78,
              left=0.150, right=0.99, top=0.865, bottom=0.185)

# ---------------- A ----------------
axA = fig.add_subplot(gs[0,0])
rows = [('Sequence',              0.805, 0.760, 0.843, SEQ),
        ('Endpoint geometry',     0.778, 0.727, 0.825, GEOM),
        ('Geometry + new static', 0.780, 0.725, 0.829, '#4FB3A6'),
        ('All new static',        0.729, 0.674, 0.779, '#7FC9BE'),
        ('Surface + electrostatics',0.713,0.658, 0.765, '#E8871E'),
        ('Contact network',       0.662, 0.604, 0.713, '#F0B27A'),
        ('Mechanical susceptibility',0.643,0.591,0.697,'#F5D0A9'),
        ('Missingness control',   0.609, 0.546, 0.670, MISS)]
y = np.arange(len(rows))[::-1]
for yi, (lab, v, l, h, c) in zip(y, rows):
    axA.barh(yi, v-0.50, left=0.50, height=0.62, color=c, edgecolor='white', lw=0.7, zorder=3)
    axA.plot([l, h], [yi, yi], color=INK, lw=0.9, zorder=5)
    axA.plot([l,l],[yi-0.15,yi+0.15], color=INK, lw=0.9, zorder=5)
    axA.plot([h,h],[yi-0.15,yi+0.15], color=INK, lw=0.9, zorder=5)
    axA.text(h+0.008, yi, f'{v:.3f}', va='center', ha='left', fontsize=6.0, color='#495057')
axA.axvline(0.5, color='#9AA0A6', lw=0.8, zorder=2)
axA.set_yticks(y); axA.set_yticklabels([r[0] for r in rows], fontsize=6.2)
axA.set_xlim(0.50, 0.92); axA.set_xticks([0.5,0.6,0.7,0.8,0.9])
axA.set_xlabel('Macro ROC AUC', labelpad=2)
axA.set_title('Representation comparison', fontsize=6.9, pad=5)
S.vgrid(axA); S.panel(axA, 'A', x=-0.72)

# ---------------- B ----------------
axB = fig.add_subplot(gs[0,1])
comp = [('Contact network',        'missingness', 0.053, -0.018, 0.130),
        ('Surface + electrostatics','missingness', 0.104,  0.031, 0.174),
        ('Mechanical susceptibility','missingness',0.034, -0.039, 0.101),
        ('All new static',         'missingness', 0.120,  0.050, 0.189),
        ('All new static',         'geometry',   -0.049, -0.104, 0.007),
        ('All new static',         'sequence',   -0.077, -0.123,-0.032),
        ('Geometry + new static',  'geometry',    0.002, -0.035, 0.037),
        ('Sequence + new static',  'sequence',   -0.054, -0.092,-0.016)]
y = np.arange(len(comp))[::-1]
for yi, (cand, base, d, l, h) in zip(y, comp):
    sig = (l > 0) or (h < 0)
    c = GEOM if d > 0 and sig else ('#C1121F' if h < 0 else '#ADB5BD')
    axB.plot([l, h], [yi, yi], color=c, lw=1.4, solid_capstyle='round', zorder=4)
    axB.scatter([d], [yi], s=32, c=c, edgecolor='white', lw=0.8, zorder=6)
axB.axvline(0, color=INK, lw=1.0, zorder=3)
axB.set_yticks(y)
axB.set_yticklabels([f'{c}\nvs {b}' for c, b, *_ in comp], fontsize=5.8, linespacing=1.12)
axB.set_xlim(-0.16, 0.22); axB.set_xlabel(r'$\Delta$ macro ROC AUC', labelpad=2)
axB.set_ylim(-0.7, len(comp)-0.05)
axB.set_title('Incremental value', fontsize=6.9, pad=5)
axB.text(0.135, len(comp)-0.45, 'adds', fontsize=5.9, color=GEOM, fontweight='bold', ha='center')
axB.text(-0.105, len(comp)-0.45, 'redundant', fontsize=5.9, color='#C1121F',
         fontweight='bold', ha='center')
S.vgrid(axB); S.panel(axB, 'B', x=-0.62)

# ---------------- C ----------------
axC = fig.add_subplot(gs[0,2])
names = ['Seq.', 'Geom.', 'Surf. +\nelec.', 'New\nstatic']
lp  = [0.812, 0.640, 0.596, 0.430]; lpl = [0.622,0.304,0.394,0.244]; lph = [0.983,0.958,0.822,0.646]
po  = [0.820, 0.644, 0.639, 0.383]; pol = [0.579,0.276,0.410,0.153]; poh = [1.000,0.944,0.944,0.595]
x = np.arange(4); w = 0.36
axC.bar(x-w/2, lp, w, color='#5C6B84', edgecolor='white', lw=0.7, zorder=3, label='Ligand-permissive (11)')
axC.bar(x+w/2, po, w, color='#B7C4D4', edgecolor='white', lw=0.7, zorder=3, label='Polymer-only (9)')
axC.errorbar(x-w/2, lp, yerr=[np.array(lp)-np.array(lpl), np.array(lph)-np.array(lp)],
             fmt='none', ecolor=INK, elinewidth=0.8, capsize=1.8, capthick=0.8, zorder=4)
axC.errorbar(x+w/2, po, yerr=[np.array(po)-np.array(pol), np.array(poh)-np.array(po)],
             fmt='none', ecolor=INK, elinewidth=0.8, capsize=1.8, capthick=0.8, zorder=4)
axC.axhline(0.5, color='#9AA0A6', lw=0.8, ls=(0,(4,2)), zorder=2)
axC.set_xticks(x); axC.set_xticklabels(names, fontsize=6.0, linespacing=1.1)
axC.set_ylim(0.10, 1.10); axC.set_ylabel('Macro ROC AUC', fontsize=6.9)
axC.set_title('Receptor-only transfer', fontsize=6.9, pad=5)
axC.legend(frameon=False, fontsize=5.7, loc='upper right', handlelength=1.1,
           handletextpad=0.4, labelspacing=0.2, borderpad=0.0)
axC.text(0.5, -0.155, 'wide intervals: current data frontier', transform=axC.transAxes,
         ha='center', fontsize=5.9, color='#495057', style='italic')
S.hgrid(axC); S.panel(axC, 'C', x=-0.34)

for ext in ('pdf','svg','png'):
    fig.savefig(f'../figures/Figure_4.{ext}', bbox_inches='tight', facecolor='white')
print('fig4 done')
