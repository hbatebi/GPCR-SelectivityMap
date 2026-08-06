import matplotlib as mpl
from matplotlib import font_manager

def apply():
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Liberation Sans', 'DejaVu Sans'],
        'font.size': 7,
        'axes.labelsize': 7.5,
        'axes.titlesize': 8,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 6.8,
        'axes.linewidth': 0.7,
        'xtick.major.width': 0.7,
        'ytick.major.width': 0.7,
        'xtick.major.size': 2.6,
        'ytick.major.size': 2.6,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'svg.fonttype': 'none',
        'figure.dpi': 200,
        'savefig.dpi': 600,
        'axes.axisbelow': True,
    })

# palette
BOUND   = '#C0396B'
FREE    = '#2E86AB'
DEFORM  = '#8D99AE'
COMB    = '#5C6B84'
SEQ     = '#6F44A8'
GEOM    = '#1B998B'
MISS    = '#B7BCC2'
ORANGE  = '#E8871E'
INK     = '#20242B'
GS, GI, GQ = '#E8871E', '#1B998B', '#6F44A8'
GRID    = '#DFE3E8'

def panel(ax, letter, x=-0.155, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9.5,
            fontweight='bold', va='bottom', ha='left', color=INK)

def vgrid(ax):
    ax.xaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)

def hgrid(ax):
    ax.yaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
