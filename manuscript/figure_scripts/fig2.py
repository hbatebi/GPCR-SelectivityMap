import sys; sys.path.insert(0,'.')
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch
import style as S; from style import *
S.apply()

fig = plt.figure(figsize=(7.0, 3.30))
gs = GridSpec(1, 2, figure=fig, width_ratios=[1.0, 1.34], wspace=0.24,
              left=0.175, right=0.99, top=0.90, bottom=0.145)

# ---------------- A ----------------
axA = fig.add_subplot(gs[0,0])
rows = [('Four-feature activation code', 0.9963, 0.9888, 1.000, ORANGE),
        ('Canonical geometry + ICL2',    0.9865, 0.9573, 1.000, '#3D5F94'),
        ('TM6 opening',                  0.9863, 0.9596, 1.000, '#5B84B1'),
        ('Metadata',                     0.8712, 0.8124, 0.9221, '#8D99AE'),
        ('ICL2 alone',                   0.5550, 0.4960, 0.6290, MISS),
        ('Feature availability',         0.5434, 0.4791, 0.6115, MISS)]
y = np.arange(len(rows))[::-1]
for yi, (lab, v, l, h, c) in zip(y, rows):
    axA.plot([l, h], [yi, yi], color=c, lw=1.6, solid_capstyle='round', zorder=3)
    axA.scatter([v], [yi], s=46, c=c, edgecolor='white', lw=1.0, zorder=5)
    axA.text(h+0.012, yi, f'{v:.3f}', va='center', ha='left', fontsize=6.4,
             fontweight='bold' if v > 0.95 else 'normal',
             color=INK if v > 0.95 else '#495057')
axA.axvline(0.5, color='#9AA0A6', lw=0.8, ls=(0,(4,2)), zorder=2)
axA.text(0.505, -0.72, 'chance', fontsize=5.9, color='#6C757D')
axA.set_yticks(y); axA.set_yticklabels([r[0] for r in rows], fontsize=6.6)
axA.set_xlim(0.44, 1.10); axA.set_xticks([0.5,0.6,0.7,0.8,0.9,1.0])
axA.set_ylim(-1.05, len(rows)-0.4)
axA.set_xlabel('ROC AUC (active vs inactive)')
axA.text(0.99, 0.055, '1,336 structures, 199 receptors',
         transform=axA.transAxes, ha='right', fontsize=6.1, color='#495057')
S.vgrid(axA); S.panel(axA, 'A', x=-0.46)

# ---------------- B ----------------
axB = fig.add_subplot(gs[0,1]); axB.set_axis_off()
img = mpimg.imread('../figures/Figure_2B_receptor_schematic.png')
axB.imshow(img, interpolation='lanczos', aspect='equal', zorder=1)
h, w = img.shape[0], img.shape[1]
axB.set_xlim(-w*0.01, w*1.01); axB.set_ylim(h*1.38, -h*0.13)

FC = ['#C0396B', '#1B998B', '#6F44A8', '#E8871E']
def arrow(n, x1, x2, y_, c):
    axB.add_patch(FancyArrowPatch((x1, y_), (x2, y_), arrowstyle='<|-|>',
                  mutation_scale=6, lw=1.3, color=c, zorder=8))
    axB.scatter([(x1+x2)/2], [y_], s=58, c=c, edgecolor='white', lw=1.0, zorder=9)
    axB.text((x1+x2)/2, y_, n, ha='center', va='center',
             fontsize=5.4 if len(n) < 3 else 4.4,
             color='white', fontweight='bold', zorder=10)

# coordinates in cropped-image space (offset x-136, y-52)
arrow('iii', 645, 840, 470, FC[2])    # PIF 5.50 - 6.44  (TM5 <-> TM6, mid-bundle)
arrow('i', 665, 915, 650, FC[0])    # TM3 3.50 - TM6 6.34 (intracellular ends)
arrow('ii', 455, 915, 778, FC[1])    # TM3 - TM7
axB.scatter([468], [640], s=58, c=FC[3], edgecolor='white', lw=1.0, zorder=9)
axB.text(468, 640, 'iv', ha='center', va='center', fontsize=5.0, color='white',
         fontweight='bold', zorder=10)

axB.text(w*0.5, -h*0.045, 'Sparse activation code', ha='center', va='bottom',
         fontsize=7.2, fontweight='bold', color=INK)
key = [('TM3 3.50 - TM6 6.34', 0.02, 1.08), ('TM3 - TM7 distance', 0.53, 1.08),
       ('PIF 5.50 - 6.44',     0.02, 1.22), ('Y7.53 displacement', 0.53, 1.22)]
for (txt, fx, fy), c, n in zip(key, [FC[0], FC[1], FC[2], FC[3]], ['i','ii','iii','iv']):
    axB.scatter([w*fx], [h*fy], s=52, c=c, edgecolor='white', lw=0.9, zorder=9)
    axB.text(w*fx, h*fy, n, ha='center', va='center',
             fontsize=5.4 if len(n) < 3 else 4.4,
             color='white', fontweight='bold', zorder=10)
    axB.text(w*fx + w*0.035, h*fy, txt, ha='left', va='center', fontsize=6.4,
             color='#343A40')
S.panel(axB, 'B', x=-0.02, y=0.99)

for ext in ('pdf','svg','png'):
    fig.savefig(f'../figures/Figure_2.{ext}', bbox_inches='tight', facecolor='white')
print('fig2 done')
