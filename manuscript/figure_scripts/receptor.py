import numpy as np
from matplotlib.patches import FancyBboxPatch, Ellipse, Circle, Polygon
from style import INK, BOUND, FREE, ORANGE

HELIX = ['#3D5F94','#4E77A6','#6592B8','#8FB3C9','#7BA6BE','#C0396B','#54789F']
NAMES = ['TM1','TM2','TM3','TM4','TM5','TM6','TM7']
YTOP, YBOT = 1.52, 0.30

#            active (bound / active-state)          inactive
G_ACT   = [(-0.74,-0.66),(-0.49,-0.41),(-0.25,-0.15),(0.00,0.05),
           (0.26,0.33),(0.50,0.66),(0.74,0.55)]
G_INACT = [(-0.74,-0.66),(-0.49,-0.41),(-0.25,-0.11),(0.00,0.05),
           (0.26,0.28),(0.50,0.40),(0.74,0.60)]

def membrane(ax, cx, scale=1.0, w=1.08, pad=0.12):
    y0, y1 = YBOT+pad, YTOP-pad
    ax.add_patch(FancyBboxPatch((cx-w*scale, y0*scale), 2*w*scale, (y1-y0)*scale,
                 boxstyle="round,pad=0,rounding_size=0.03", fc='#F6F2E7',
                 ec='none', zorder=0))
    for yy in (y0, y1):
        ax.plot([cx-w*scale, cx+w*scale], [yy*scale, yy*scale],
                color='#DFD5B8', lw=1.0, zorder=0.5)

def bundle(ax, cx, active=True, scale=1.0, lw=5.4, label_tm=(), zorder=3):
    g = G_ACT if active else G_INACT
    for i, (xt, xb) in enumerate(g):
        ax.plot([cx+xt*scale, cx+xb*scale], [YTOP*scale, YBOT*scale],
                color=HELIX[i], lw=lw, solid_capstyle='round', zorder=zorder)
        if NAMES[i] in label_tm:
            ax.text(cx+xb*scale+0.10*scale, YBOT*scale-0.02*scale, NAMES[i],
                    ha='left', va='center', fontsize=5.8,
                    color=HELIX[i], fontweight='bold', zorder=9)
    ax.plot([cx+0.60*scale, cx+0.95*scale], [(YBOT-0.06)*scale, (YBOT-0.09)*scale],
            color='#54789F', lw=lw*0.8, solid_capstyle='round', zorder=zorder)

def ligand(ax, cx, scale=1.0, color=ORANGE, label=False):
    th = np.linspace(0, 2*np.pi, 7, endpoint=False)+0.4
    r = 0.085*scale
    ax.add_patch(Polygon(np.c_[cx+r*np.cos(th), (YTOP-0.30)*scale+r*np.sin(th)],
                 closed=True, fc=color, ec='white', lw=0.8, zorder=6))
    if label:
        ax.text(cx-0.16*scale, (YTOP-0.30)*scale, 'ligand', ha='right', va='center',
                fontsize=5.8, color='#9A6316')

def gprotein(ax, cx, scale=1.0):
    ax.plot([cx+0.10*scale, cx+0.02*scale], [-0.20*scale, (YBOT+0.16)*scale],
            color='#8A2846', lw=5.2, solid_capstyle='round', zorder=2)
    ax.text(cx+0.19*scale, 0.10*scale, r'$\alpha$5', fontsize=5.8, color='#8A2846',
            fontweight='bold', ha='left', va='center', zorder=9)
    ax.add_patch(Ellipse((cx+0.06*scale, -0.52*scale), 1.02*scale, 0.50*scale,
                 fc=BOUND, ec='white', lw=1.0, zorder=4))
    ax.text(cx+0.06*scale, -0.52*scale, r'G$\alpha$', ha='center', va='center',
            fontsize=6.6, color='white', fontweight='bold', zorder=5)
    ax.add_patch(Circle((cx-0.68*scale, -0.60*scale), 0.24*scale,
                 fc='#E08AA6', ec='white', lw=1.0, zorder=4))
    ax.text(cx-0.68*scale, -0.60*scale, r'G$\beta$', ha='center', va='center',
            fontsize=5.6, color='white', fontweight='bold', zorder=5)
    ax.add_patch(Circle((cx-0.90*scale, -0.30*scale), 0.155*scale,
                 fc='#F2C2D1', ec='white', lw=1.0, zorder=4))
    ax.text(cx-0.90*scale, -0.30*scale, r'G$\gamma$', ha='center', va='center',
            fontsize=4.8, color='#8A2846', fontweight='bold', zorder=5)
